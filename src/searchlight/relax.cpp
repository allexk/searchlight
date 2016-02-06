/*
 * relax.cpp
 *
 *  Created on: Jan 3, 2016
 *      Author: akalinin
 */

#include "relax.h"
#include "searchlight.h"

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
    log4cxx::Logger::getLogger("searchlight.relax"));

// Relaxable constraint tags defined
const char RelaxableConstraint::BetweenConstTag[] = "RelaxableBetween";
const char RelaxableConstraint::LessEqConstTag[] = "RelaxableLEQ";
const char RelaxableConstraint::GreaterEqConstTag[] = "RelaxableGEQ";

// Output to a stream
std::ostream &operator<<(std::ostream &os, const Relaxator::FailReplay &fr) {
    os << "Fail replay:";
    for (const auto &fc: fr.failed_const_) {
        os << ' ' << fc.rl_ << "<=" << fc.const_id_ << "<=" << fc.rh_;
    }
    os << ": RD=" << fr.relax_degree_;
    return os;
}

void Relaxator::RegisterConstraint(const std::string &name, size_t solver_id,
		RelaxableConstraint *constr, int64 max_l, int64 max_h) {
	assert(constr && solver_id < solvers_num_);

	// Find the constraint info
	if (const_to_id_.find(name) == const_to_id_.end()) {
		// Interval values
		int64 l, h;
		constr->GetParams(l, h);

		// Determine the constraint type
		IntervalConstraint::ConstraintType constr_type;
		if (h == kint64max) {
			constr_type = IntervalConstraint::ConstraintType::GEQ;
			max_h = kint64max;
		} else if (l == kint64min) {
			constr_type = IntervalConstraint::ConstraintType::LEQ;
			max_l = kint64min;
		} else {
			constr_type = IntervalConstraint::ConstraintType::BET;
		}
		/*
		 *  Adjust max_l and max_h to avoid negative distances. In this case,
		 *  however, the relaxation along the corresponding bound won't be
		 *  possible.
		 */
		if (max_l > l) {
			max_l = l;
		}
		if (max_h < h) {
			max_h = h;
		}

		// Create a new constraint info
		orig_consts_.emplace_back(IntervalConstraint(constr_type, l, h),
				solvers_num_, max_l, max_h);
		const_to_id_[name] = orig_consts_.size() - 1;
	}
	// Id now exists
	const size_t const_id = const_to_id_[name];
	constr->SetId(const_id);
    ConstraintInfo &const_info = orig_consts_[const_id];
	const_info.solver_const_[solver_id] = constr;
	LOG4CXX_INFO(logger, " Registered constraint (" << name << ") for solver="
	    << solver_id <<" with relaxation interval (" << max_l << ','
	    << max_h << "}, id=" << const_id);
}

int64 Relaxator::ComputeNewReplayInterval(
		FailReplay::FailedConstraint &failed_const, double max_relax_dist,
		int rel_pos) {
	const ConstraintInfo &ci = orig_consts_[failed_const.const_id_];
	const int64 max_relax_dist_int =
			ci.MaxRelaxDist(rel_pos == -1)  * max_relax_dist;

	// Compute the distance
	const int64 min_r_dist = failed_const.rl_;
	const int64 max_r_dist = failed_const.rh_;
	int64 r_dist;
	if (max_relax_dist_int > max_r_dist) {
		r_dist = max_r_dist;
	} else if (max_relax_dist_int >= min_r_dist) {
		r_dist = max_relax_dist_int;
	} else {
		// Cannot relax
		return -1;
	}
	// Compute new interval
	if (rel_pos == 1) {
		failed_const.rl_ = ci.int_.Min();
		failed_const.rh_ = ci.int_.Max() + r_dist;
	} else {
		failed_const.rl_ = ci.int_.Min() - r_dist;
		failed_const.rh_ = ci.int_.Max();
	}
	return r_dist;
}

void Relaxator::RegisterFail(size_t solver_id) {
	// Compute relaxation degree: first pass, determine violated constraints
	std::vector<int> saved_rel_pos;
	FailReplay replay;
	for (size_t i = 0; i < orig_consts_.size(); ++i) {
		const ConstraintInfo &ci = orig_consts_[i];
		const IntExpr *c_expr = ci.solver_const_[solver_id]->GetExpr();
		const int64 cmin = c_expr->Min();
		const int64 cmax = c_expr->Max();
		// if cmin > cmax, that means the constraint cannot be fulfilled ever
		if (cmin > cmax) {
			return;
		}
		// Check if we're violating the constraint
		int64 min_d, max_d;
		const int rel_pos = ci.int_.MinMaxDistance(cmin, cmax, min_d, max_d);
		if (rel_pos != 0) {
			// Violation
			saved_rel_pos.push_back(rel_pos);
			replay.failed_const_.push_back({i, cmin, cmax, min_d, max_d});
		}
	}
	if (saved_rel_pos.empty()) {
	    // No violations: might be a "custom" fail
	    return;
	}
	// Second pass: compute the degree
	double relax_dist = 0.0;
	const size_t viol_const_num = replay.failed_const_.size();
	const double max_relax_dist = MaxUnitRelaxDistance(viol_const_num);
	for (size_t i = 0; i < replay.failed_const_.size(); ++i) {
		FailReplay::FailedConstraint &failed_const = replay.failed_const_[i];
		const int64 r_dist = ComputeNewReplayInterval(failed_const,
				max_relax_dist, saved_rel_pos[i]);
		if (r_dist == -1) {
			// Cannot relax: ignore the fail
			return;
		}
		const ConstraintInfo &ci = orig_consts_[failed_const.const_id_];
		relax_dist = std::max(relax_dist,
				double(r_dist) / ci.MaxRelaxDist(saved_rel_pos[i] == -1));
	}
	replay.relax_degree_ = RelaxDistance(relax_dist, viol_const_num);
	// Save the variables for the future replay
	replay.saved_vars_ = sl_.GetSLSolver(solver_id).GetCurrentVarDomains();
	// Debug print
	LOG4CXX_DEBUG(logger, "New replay registered: " << replay);
	// Put the replay into the queue
	{
		std::lock_guard<std::mutex> lock{mtx_};
		fail_replays_.push(std::move(replay));
	}
}

double Relaxator::ComputeResultRelaxationDegree(
        const Int64Vector &rel_result) const {
    assert(orig_consts_.size() == rel_result.size());
    double max_relax = 0;
    size_t viol_consts = 0;
    for (size_t i = 0; i < rel_result.size(); ++i) {
        const double rel_relax = orig_consts_[i].RelRelaxDist(rel_result[i]);
        if (rel_relax != 0.0) {
            // Comparison with 0 is fine here
            max_relax = std::max(rel_relax, max_relax);
            viol_consts++;
        }
    }
    return RelaxDistance(max_relax, viol_consts);
}

void Relaxator::ReportResult(double rd) {
    std::lock_guard<std::mutex> lock{mtx_};
    if (top_results_.size() < res_num_) {
        top_results_.push(rd);
    } else {
        double min_rd = top_results_.top();
        if (min_rd >= rd) {
            top_results_.pop();
            top_results_.push(rd);
            min_rd = top_results_.top();
            if (min_rd + 0.001 < lrd_.load(std::memory_order_relaxed)) {
                // change LRD (use some EPSILON to avoid near close changes)
                lrd_.store(min_rd, std::memory_order_relaxed);
                LOG4CXX_DEBUG(logger, "LRD changed to " << min_rd);
            }
        }
    }
}
} /* namespace searchlight */

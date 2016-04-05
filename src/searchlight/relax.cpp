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

// Model tags
const char RelaxableConstraint::ModelIDTag[] = "constraint_id";


// Output to a stream
std::ostream &operator<<(std::ostream &os, const Relaxator::FailReplay &fr) {
    os << "Fail replay:";
    for (const auto &fc: fr.failed_const_) {
        os << ' ' << fc.rl_ << "<=" << fc.const_id_ << "<=" << fc.rh_;
    }
    os << ": best_RD=" << fr.best_relax_degree_ << " worst_RD="
            << fr.worst_relax_degree_;
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
		/*
		 * Check if the constraint can be relaxed at all. If not, ignore.
		 */
		if (l == max_l && h == max_h) {
		    LOG4CXX_INFO(logger, "Constraint (" << name <<
		            ") cannot be relaxed due to user's values. Ignoring...");
		    return;
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

bool Relaxator::ComputeNewReplayInterval(FailReplay &replay,
		FailReplay::FailedConstraint &failed_const,
		double max_relax_dist) const {
	const ConstraintInfo &ci = orig_consts_[failed_const.const_id_];
	// Maximum relaxation distance for the constraint
	const int64 max_absolute_relax_dist_int =
	        ci.MaxRelaxDist(failed_const.rel_pos_ == -1);
	// Maximum relaxation distance with respect to the current lrd
	const int64 max_relax_dist_int =
	        max_absolute_relax_dist_int  * max_relax_dist;

    /*
     *  Compute new interval. We have a choice here. We can just extend the
     *  original interval or create a new one, based on the current (fail)
     *  values. Instead, we store the minimum and maximum
     *  relaxation distance, and the interval can be restored later in the way
     *  we wish.
     */
	if (max_relax_dist_int >= failed_const.rl_ &&
	        max_relax_dist_int < failed_const.rh_) {
	    // Tighten the maximum relaxation
	    failed_const.rh_ = max_relax_dist_int;
	} else if (max_relax_dist_int < failed_const.rl_) {
		// Cannot relax with respect to this fail
		return false;
	}
	replay.best_relax_degree_ = std::max(replay.best_relax_degree_,
	    double(failed_const.rl_) / max_absolute_relax_dist_int);
    replay.worst_relax_degree_ = std::max(replay.worst_relax_degree_,
        double(failed_const.rh_) / max_absolute_relax_dist_int);
    return true;
}

void Relaxator::RegisterFail(size_t solver_id) {
	// Compute relaxation degree: first pass, determine violated constraints
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
			replay.failed_const_.push_back({i, rel_pos, min_d, max_d});
		}
	}
	if (replay.failed_const_.empty()) {
	    // No violations: might be a "custom" fail
	    return;
	}
	// Second pass: compute the degree
	if (!ComputeFailReplayRelaxation(replay)) {
        // Cannot relax: ignore the fail
	    return;
	}
	// Normalize relax degrees
    const size_t viol_const_num = replay.failed_const_.size();
	replay.best_relax_degree_ = RelaxDistance(replay.best_relax_degree_,
	        viol_const_num);
    replay.worst_relax_degree_ = RelaxDistance(replay.worst_relax_degree_,
            viol_const_num);
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

double Relaxator::ComputeViolSpecBestRelaxationDegree(
    const Int64Vector &vc) const {

    double max_relax = 0;
    size_t violated = 0;
    for (size_t i = 0; i < vc.size(); i += 3) {
        const int64 cid = vc[i];
        const int64 l = vc[i + 1];
        const int64 h = vc[i + 2];
        const double rd = orig_consts_[cid].MinRelRelaxDist(l, h);
        if (rd != 0.0) {
            max_relax = std::max(max_relax, rd);
            ++violated;
        }
    }
    return RelaxDistance(max_relax, violated);
}

void Relaxator::ReportResult(double rd) {
    std::lock_guard<std::mutex> lock{mtx_};
    bool lrd_might_change = false;
    if (top_results_.size() < res_num_) {
        top_results_.push(rd);
        lrd_might_change = top_results_.size() == res_num_;
    } else if (top_results_.top() > rd) {
        top_results_.pop();
        top_results_.push(rd);
        lrd_might_change = true;
    }
    if (lrd_might_change) {
        const double min_rd = top_results_.top();
        if (min_rd + 0.001 < lrd_.load(std::memory_order_relaxed)) {
            // change LRD (use some EPSILON to avoid near close changes)
            lrd_.store(min_rd, std::memory_order_relaxed);
            LOG4CXX_DEBUG(logger, "LRD changed to " << min_rd);
        }
    }
}

bool Relaxator::ComputeFailReplayRelaxation(FailReplay &replay) const {
    const size_t viol_const_num = replay.failed_const_.size();
    assert(viol_const_num);
    // Maximum relative relaxation based on the current LRD
    const double max_relax_dist = MaxUnitRelaxDistance(viol_const_num);
    if (max_relax_dist <= 0.0) {
        // Cannot relax with this number of violations
        return false;
    }
    replay.best_relax_degree_ = replay.worst_relax_degree_ = 0;
    for (auto &failed_const: replay.failed_const_) {
        if (!ComputeNewReplayInterval(replay, failed_const, max_relax_dist)) {
            // Cannot relax: ignore the fail
            return false;
        }
    }
    return true;
}

Int64Vector Relaxator::ViolConstSpec(const FailReplay &replay) const {
    Int64Vector res;
    res.reserve(replay.failed_const_.size() * 3);
    for (const auto &fc: replay.failed_const_) {
        const IntervalConstraint &ic = orig_consts_[fc.const_id_].int_;
        res.push_back(fc.const_id_);
        int64 l, h;
        if (fc.rel_pos_ == 1) {
            l = ic.Max() + fc.rl_;
            h = ic.Max() + fc.rh_;
        } else {
            l = ic.Min() - fc.rh_;
            h = ic.Min() - fc.rl_;
        }
        res.push_back(l);
        res.push_back(h);
    }
    return res;
}

Int64Vector Relaxator::GetMaximumRelaxationVCSpec() const {
    Int64Vector res(orig_consts_.size() * 3);
    // Maximum distance is reached when only 1 constraint is violated
    const double max_relax_dist = MaxUnitRelaxDistance(1);
    if (max_relax_dist <= 0.0) {
        // No violations possible
        res.clear();
        return res;
    }
    for (size_t i = 0; i < orig_consts_.size(); ++i) {
        res[3 * i] = i;
        orig_consts_[i].MaxRelaxDist(max_relax_dist, res[3 * i + 1],
                res[3 * i + 2]);
    }
    return res;
}

FailCollectorMonitor::FailCollectorMonitor(Solver &solver,
        const SearchlightSolver &sl_solver,
        Relaxator &rel) :
    SearchMonitor(&solver),
    relaxator_(rel),
    sl_solver_(sl_solver),
    solver_id_(sl_solver.GetLocalID()) {}

void FailCollectorMonitor::BeginFail() {
    if (!sl_solver_.LastFailCustom()) {
        relaxator_.RegisterFail(solver_id_);
    }
}

} /* namespace searchlight */

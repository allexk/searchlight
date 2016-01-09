/*
 * relax.cpp
 *
 *  Created on: Jan 3, 2016
 *      Author: akalinin
 */

#include "relax.h"

namespace searchlight {

// Relaxable constraint tags defined
const char RelaxableConstraint::BetweenConstTag[] = "RelaxableBetween";
const char RelaxableConstraint::LessEqConstTag[] = "RelaxableLEQ";
const char RelaxableConstraint::GreaterEqConstTag[] = "RelaxableGEQ";

void Relaxator::RegisterConstraint(const std::string &name, size_t solver_id,
		RelaxableConstraint *constr, IntExpr *expr, int64 max_l, int64 max_h) {
	assert(constr && expr && solver_id < solvers_num_);

	// Find the constraint info
	const auto it = const_to_id_.find(name);
	ConstraintInfo *const_info;
	if (it == const_to_id_.end()) {
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
		// Adjust max_l and max_h to allow at least the smallest relaxation
		if (max_l > l) {
			max_l = l - 1;
		}
		if (max_h < h) {
			max_h = h + 1;
		}

		// Create a new constraint info
		orig_consts_.emplace_back(IntervalConstraint(constr_type, l, h),
				solvers_num_, max_l, max_h);
		const_to_id_[name] = orig_consts_.size() - 1;
		const_info = &orig_consts_.back();
	} else {
		const_info = &orig_consts_[it->second];
	}
	const_info->solver_const_[solver_id] = constr;
	const_info->solver_exprs_[solver_id] = expr;
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
		const IntExpr *c_expr = ci.solver_exprs_[solver_id];
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
	// Second pass: compute the degree
	assert(!saved_rel_pos.empty());
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
	// Put the replay into the queue
	{
		std::lock_guard<std::mutex> lock{mtx_};
		fail_replays_.push(std::move(replay));
	}
}
} /* namespace searchlight */

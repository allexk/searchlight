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

std::ostream &operator<<(std::ostream &os, const Relaxator::RelaxatorStats &rs) {
    os << "Relaxator stats:\n";
    const float total_secs =
            std::chrono::duration_cast<std::chrono::duration<double>>(
                   rs.total_fail_time_).count();
    const float total_success_secs =
            std::chrono::duration_cast<std::chrono::duration<double>>(
                   rs.total_success_fail_time_).count();
    os << "\tTotal fail catch time: " << total_secs << "s\n";
    os << "\tTotal success fail catch time: " << total_success_secs << "s\n";
    os << "\tTotal fails caught: " <<
            rs.total_fails_caught_.load(std::memory_order_relaxed) << '\n';
    os << "\tTotal fails registered: " << rs.total_fails_registered_ << '\n';
    os << "\tTotal fails heuristic-retried: " <<
            rs.total_fails_heur_retried_ << '\n';
    os << "\tTotal constraints re-guessed: " <<
            rs.total_const_reguessed_.load(std::memory_order_relaxed) << '\n';
    os << "\tTotal fails replayed: " << rs.total_fails_replayed_ << '\n';
    return os;
}

Relaxator::Relaxator(Searchlight &sl, size_t solvers,
    double dist_w, size_t res_num) :

    sl_(sl),
    solvers_num_(solvers),
    solver_info_(solvers),
    distance_weight_(dist_w),
    res_num_(res_num),
    register_heur_(RegisterHeuristic::ALL),
    replay_relax_(ReplayRelaxation::VIOLATED),
    sort_method_(ReplaySortMethod::BEST) {

    // Default fail heuristic
    const std::string heur = sl.GetConfig().get("relax.heur", "all");
    if (heur == "guess") {
        register_heur_ = RegisterHeuristic::GUESS;
    } else if (heur == "guess-all") {
        register_heur_ = RegisterHeuristic::GUESS_ALL;
    } else {
        LOG4CXX_ERROR(logger,
                "Unknown register heuristic, defaulting to ALL...");
    }
    // Replay relaxation
    const std::string rr = sl.GetConfig().get("relax.replay", "viol");
    if (rr == "all") {
        replay_relax_ = ReplayRelaxation::ALL;
    }
    // Replay sorting
    const std::string rep_sort = sl.GetConfig().get("relax.sort", "best");
    if (rep_sort == "worst") {
        sort_method_ = ReplaySortMethod::WORST;
    } else if (rep_sort == "prod") {
        sort_method_ = ReplaySortMethod::PROD;
    } else if (rep_sort == "time") {
        sort_method_ = ReplaySortMethod::TIME;
    } else {
        LOG4CXX_ERROR(logger,
                "Unknown replay sorting method, defaulting to BEST...");
    }
    fail_replays_ = decltype(fail_replays_)(ReplaySort(sort_method_), {});

    // Do we save UDFS?
    save_udfs_for_replay_ = sl.GetConfig().get("relax.save_udfs", false);

    // Replay relaxation degree
    replay_rd_ = sl.GetConfig().get("relax.replay_rd", 1.0);
    if (replay_rd_ > 1.0) {
        replay_rd_ = 1.0;
    }

    // Initial LRD
    lrd_ = sl.GetConfig().get("relax.lrd", 1.0);
    if (lrd_ != 1.0) {
        LOG4CXX_WARN(logger, "Initial LRD set to " << lrd_.load());
    }
}

Relaxator::~Relaxator() {
    std::ostringstream os;
    os << "Pre-replay stage stats:\n";
    os << stats_[0];
    os << "Replay stage stats:\n";
    os << stats_[1];
    LOG4CXX_INFO(logger, os.str());
}

bool Relaxator::GetFailReplay(size_t solver_id,
        std::vector<IntVarDomainInfo> &dom_info, Int64Vector &vc_spec,
        UDFStates &udfs) {
    std::lock_guard<std::mutex> lock{mtx_};
    const double lrd = lrd_.load(std::memory_order_relaxed);
    while (!fail_replays_.empty()) {
        // Const cast is okay here; we get rid of the top object right away
        FailReplay &fr = const_cast<FailReplay &>(fail_replays_.top());
        if (fr.best_relax_degree_ <= lrd) {
            if (fr.worst_relax_degree_ > lrd) {
                // Tighten the relaxation if possible
                ComputeFailReplayRelaxation(fr);
            }
            vc_spec = ViolConstSpec(fr);
            // Get the info from fr; can move here since it'll be destroyed
            dom_info = std::move(fr.saved_vars_);
            udfs = std::move(fr.saved_udfs_);
            // For the last replay we need only info about failed constraints
            SolverReplayInfo &solver_info = solver_info_[solver_id];
            if (register_heur_ == RegisterHeuristic::GUESS) {
                for (auto &fc: fr.failed_const_) {
                    solver_info.replay_.failed_consts_.emplace(fc.const_id_,
                            std::move(fc));
                }
            }
            stats_[1].total_fails_replayed_++;
            solver_info.in_replay_ = true;
            fail_replays_.pop();
            return true;
        } else {
            // Cannot replay due to LRD, ignore...
            fail_replays_.pop();
            if (sort_method_ == ReplaySortMethod::BEST) {
                // Actually can clear the queue, since all others > lrd
                fail_replays_ = decltype(fail_replays_)();
            }
        }
    }
    return false;
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
		 * Note, we have to register the constraint even if it cannot be
		 * relaxed, since at RegisterFail we have to understand if this
		 * constraint is satisfied or not to avoid replays that cannot be
		 * remedied.
		 *
		 * TODO: this is a larger issue of mixing relaxable and usual
		 * constraints in the same query. If the "common" and relaxable one
		 * fail at the same time, the fail might get registered, but
		 * will fail again immediately, resulting in an infinite loop.
		 * Should try to avoid this. For now the only way to do this is to
		 * make all constraints Relaxable. Another way is to record the
		 * last replay for every solver and check if the new fail is equal
		 * to the last replay -- if the replay was due to the non-relaxable
		 * constraint, it will fail immediately at the replay and the new
		 * replay will equal to the last one.
		 */
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

void Relaxator::UpdateTimeStats(RelaxatorStats &stats,
        const std::chrono::steady_clock::time_point &start_time,
        bool success) {
    // Stop time
    const auto end_time = std::chrono::steady_clock::now();
    stats.total_fail_time_ +=
            std::chrono::duration_cast<std::chrono::microseconds>(
                    end_time - start_time);
    if (success) {
        stats.total_success_fail_time_ +=
                std::chrono::duration_cast<std::chrono::microseconds>(
                        end_time - start_time);
    }
}

void Relaxator::RegisterFail(size_t solver_id, RegisterHeuristic rh) {
	// Compute relaxation degree: first pass, determine violated constraints
    const SolverReplayInfo &solver_info = solver_info_[solver_id];
    RelaxatorStats &rel_stats = stats_[solver_info.GetStatsFrame()];
    rel_stats.total_fails_caught_.fetch_add(1, std::memory_order_relaxed);
    // starting the timer
    const auto reg_start_time = std::chrono::steady_clock::now();
	FailReplay replay;
	if (rh == RegisterHeuristic::GUESS) {
	    // With GUESS heuristic we don't want to recompute the functions
	    sl_.GetSLSolver(solver_id).SetAdapterMode(Adapter::DUMB, true);
	}
	bool cannot_relax = false;
	for (size_t i = 0; i < orig_consts_.size(); ++i) {
		const ConstraintInfo &ci = orig_consts_[i];
		const IntExpr *c_expr = ci.solver_const_[solver_id]->GetExpr();
		const int64 cmin = c_expr->Min();
		const int64 cmax = c_expr->Max();
		// if cmin > cmax, that means the constraint cannot be fulfilled ever
		if (cmin > cmax) {
		    cannot_relax = true;
		    break;
		}
		// Check if we're violating the constraint
		int64 min_d, max_d;
		const int rel_pos = ci.int_.MinMaxDistance(cmin, cmax, min_d, max_d);
		// Check if we can relax
		if (rel_pos != 0) {
	        // Check if we can relax
		    if (!ci.CanRelax(rel_pos < 0, min_d)) {
		        // This is a shortcut to avoid checking all constraints
	            cannot_relax = true;
	            break;
		    }
			// Violation
			replay.failed_const_.push_back({i, rel_pos, min_d, max_d});
		} else if (register_heur_ == RegisterHeuristic::GUESS &&
		        solver_info.in_replay_) {
		    // With replay-fails on GUESS we might have more info
		    const auto it = solver_info.replay_.failed_consts_.find(i);
		    if (it != solver_info.replay_.failed_consts_.end()) {
		        rel_stats.total_const_reguessed_.fetch_add(
		                1, std::memory_order_relaxed);
		        replay.failed_const_.push_back(it->second);
		    }
		}
	}
    if (rh == RegisterHeuristic::GUESS) {
        sl_.GetSLSolver(solver_id).RestoreAdapterMode();
    }
	if (replay.failed_const_.empty() || cannot_relax) {
	    // No violations detected or cannot relax
        UpdateTimeStats(rel_stats, reg_start_time, false);
        if (rh == RegisterHeuristic::GUESS && !cannot_relax) {
            // Try with the full heuristic to avoid losing fails
            rel_stats.total_fails_caught_.fetch_sub(1, std::memory_order_relaxed);
            rel_stats.total_fails_heur_retried_.fetch_add(1,
                    std::memory_order_relaxed);
            LOG4CXX_DEBUG(logger, "Retrying the fail with the ALL heuristic");
            RegisterFail(solver_id, RegisterHeuristic::ALL);
        }
	    return;
	}
	// Second pass: compute the degree
	if (!ComputeFailReplayRelaxation(replay)) {
        // Cannot relax: ignore the fail
        UpdateTimeStats(rel_stats, reg_start_time, false);
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
	// Save the UDFs
	if (save_udfs_for_replay_) {
	    sl_.GetSLSolver(solver_id).SaveUDFs(replay.saved_udfs_);
	}
	// Debug print
	LOG4CXX_DEBUG(logger, "New replay registered: " << replay);
	// Put the replay into the queue
	{
		std::lock_guard<std::mutex> lock{mtx_};
		replay.timestamp_ = fail_stamp_++;
		fail_replays_.push(std::move(replay));
		rel_stats.total_fails_registered_++;
        UpdateTimeStats(rel_stats, reg_start_time, true);
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
    const bool relaxing_all = replay_relax_ == ReplayRelaxation::ALL;
    const size_t violed_const_num = replay.failed_const_.size();
    const size_t relaxing_const_num = relaxing_all ?
            orig_consts_.size() : violed_const_num;
    std::vector<bool> relaxed_const_bs(orig_consts_.size());
    res.reserve(relaxing_const_num * 3); // triplets: id, left, right
    for (const auto &fc: replay.failed_const_) {
        const IntervalConstraint &ic = orig_consts_[fc.const_id_].int_;
        relaxed_const_bs[fc.const_id_] = true;
        res.push_back(fc.const_id_);
        int64 l, h;
        // Right distance might be modified by the replay RD.
        const int64 right_dist = fc.rl_ + (fc.rh_ - fc.rl_) * replay_rd_;
        if (fc.rel_pos_ == 1) {
            l = ic.Max() + fc.rl_;
            h = ic.Max() + right_dist;
        } else {
            l = ic.Min() - right_dist;
            h = ic.Min() - fc.rl_;
        }
        res.push_back(l);
        res.push_back(h);
    }
    // If we use ALL relaxation, relax the rest
    if (relaxing_all && relaxing_const_num > violed_const_num) {
        /*
         *  The distance is computed based on the violated constraints! It
         *  might be larger than necessary, but we avoid re-replays this
         *  way: such relaxations cannot be relaxed any further.
         *
         *  We can safely use +1 here, since in case such a constraint really
         *  fails, it's going to be in addition to others.
         *
         *  We also take replay_rd into consideration here, since it
         *  governs how we want to relax replays.
         */
        const double max_relax_dist =
                MaxUnitRelaxDistance(violed_const_num + 1) * replay_rd_;
        // linear search over the bitset is fine here
        for (size_t i = 0; i < relaxed_const_bs.size(); ++i) {
            if (!relaxed_const_bs[i]) {
                int64_t l, h;
                orig_consts_[i].MaxRelaxDist(max_relax_dist, l, h);
                res.push_back(i);
                res.push_back(l);
                res.push_back(h);
            }
        }
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

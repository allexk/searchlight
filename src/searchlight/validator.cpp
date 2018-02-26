/*
 * Copyright 2014, Brown University, Providence, RI.
 *
 *                         All Rights Reserved
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose other than its incorporation into a
 * commercial product is hereby granted without fee, provided that the
 * above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation, and that the name of Brown University not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific, written prior permission.
 *
 * BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
 * INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
 * PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
 * ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

/**
 * @file validator.cpp
 * The implementation of the validator.
 *
 * @author Alexander Kalinin
 */

#include "searchlight.h"
#include "validator.h"
#include "searchlight_task.h"

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.validator"));

// Output validator stats to a stream
std::ostream &operator<<(std::ostream &os, const Validator::Stats &vs) {
    os << "Validator stats:\n";
    os << "\tTotal canidates: " <<
            vs.total_local_cands_ + vs.total_remote_cands_ <<
            " (local: " << vs.total_local_cands_ <<
            ", remote: " << vs.total_remote_cands_ << ")\n";
    os << "\tForwarded: " << vs.total_forw_ << '\n';
    os << "\tRe-forwarded: " << vs.total_reforw_ << '\n';
    os << "\tLocal checks: " << vs.local_check_ << '\n';
    os << "\tRemote checks: " << vs.remote_check_ << '\n';
    os << "\tRelax pre-filter: " << vs.relax_pre_filtered_ << '\n';
    os << "\tRelax post-filter: " << vs.relax_post_filtered_ << '\n';
    return os;
}

/**
 * A decision that just restores the given assignment causing the corresponding
 * variables to take the values from the assignment. In its right branch it
 * does nothing, thus creating a right-deep tree.
 */
class RestoreAssignment : public Decision {
public:
    /**
     * Constructs a restore assignment decision.
     *
     * @param asgn the assignment to restore
     * @param relaxed_constrs relaxed constraints
     * @param rcs relaxable constraints
     */
    RestoreAssignment(Assignment *asgn, const Int64Vector &relaxed_constrs,
                      const RelaxableConstraints &rcs) :
        asgn_(asgn),
        relaxed_constrs_(relaxed_constrs),
        constrs_(rcs) {}

    /**
     * The left branch. In this case it restores the assignment.
     *
     * @param s the solver
     */
    virtual void Apply(Solver* const s) override {
        // Set relaxation constraints
        for (size_t i = 0; i < relaxed_constrs_.size(); i += 3) {
            const int64_t cid = relaxed_constrs_[i];
            const int64_t l_param = relaxed_constrs_[i + 1];
            const int64_t h_param = relaxed_constrs_[i + 2];
            constrs_[cid]->SetNewParams(l_param, h_param);
        }
        // Check the assignment
        asgn_->Restore();
    }

    /**
     * The right branch. It does nothing.
     *
     * @param s the solver
     */
    virtual void Refute(Solver* const s) override {}

    /**
     * Returns a string for debug logging.
     *
     * @return a debug string
     */
    virtual std::string DebugString() const override {
        return "RestoreAssignmentDecision";
    }

private:
    // The assignment to restore
    Assignment *asgn_;

    // Relaxed constraints
    const Int64Vector &relaxed_constrs_;

    // Relaxable constraints
    const RelaxableConstraints &constrs_;
};

/**
 * This is a DB that just generates decisions to restore the current assignment,
 * causing the corresponding variables to take the assignment's values.
 *
 * The main idea is that the validator is going to run this DB in a thread
 * until it is stopped. Then it will produce NULL, which will fail the search.
 * For this, we generate a right-deep tree: left node checks an assignment,
 * the right one does nothing.
 *
 * Since we have to produce a leaf after every successful restoration, this
 * db produces alternate decisions: after a successful restore a NULL is
 * generated that cuts the left branch and signals a leaf. This is mostly
 * taken care of by the solver itself; we just need to return NULL.
 *
 * Note that we also inherit from SearchMonitor since we need its functionality.
 * First of all, we need to restart the search to perform garbage collecting.
 * Secondly, we need to detect when the current assignment was validated
 * successfully to produce NULL and, thus, signal a leaf to the solver. This
 * is done via AfterDecision().
 */
class RestoreAssignmentBuilder : public DecisionBuilder {
private:
    /*
     *  Contains information about an action the DB issued.
     *
     *  The definition was put here, since Eclipse's codan freaks out otherwise.
     */
    struct Action {
        // Type of action
        enum class Type {
            NOP, // Nothing: no-operation
            SIM, // Simulation
            LOCAL_CHECK, // Checking local candidate
            REMOTE_CHECK // Checking remote candidate
        };

        // Type of action
        Type type_ = Type::NOP;

        // Forward id for remote solutions
        int64_t forward_id_;
    };
public:

    /**
     * Creates a new restore assignment builder.
     *
     * @param validator servicing validator
     * @param adapter adapter used for access
     * @param s the solver
     * @param rcs relaxable constraints
     * @param restart_period validator restart period
     * @param prototype prototype assignment to generate validations
     * @param track_prototype tracking variables prototype
     * @param relax_prototype relaxable constraints expressions prototype
     * @param asgns initial assignments to validate
     * @param slave true, if this is a slave builder for a helper
     */
    explicit RestoreAssignmentBuilder(Validator &validator,
            const AdapterPtr &adapter, Solver *s,
            const RelaxableConstraints &rcs,
            int restart_period, const Assignment *prototype,
			const Assignment &track_prototype,
			const Assignment &relax_prototype,
            CandidateVector &&asgns,
            bool slave)
        : validator_(validator),
          adapter_{adapter},
          constrs_(rcs),
          asgns_{std::move(asgns)},
          action_succeeded_(false),
          aux_monitor_(*this, s, restart_period),
          last_asgn_{new Assignment{prototype}},
		  track_vars_asgn_{&track_prototype},
		  relax_asgn_{&relax_prototype},
          slave_{slave} {}


    /**
     * Establishes additional monitors for the search.
     *
     * @param solver the solver
     * @param extras vector to add monitors to
     */
    virtual void AppendMonitors(Solver* const solver,
                                std::vector<SearchMonitor*>* const extras) {
        extras->push_back(&aux_monitor_);
    }

    /**
     * Produces a decision to assign values to variables, which effectively
     * validates the assignment.
     *
     * @param solver the current solver
     * @return the next decision or nullptr to signify a leaf
     */
    virtual Decision* Next(Solver* const solver) {
        if (last_action_.type_ != Action::Type::NOP) {
            // We performed some action: should handle the result
            if (last_action_.type_ == Action::Type::LOCAL_CHECK) {
                validator_.stats_.local_check_.fetch_add(
                        1, std::memory_order_relaxed);
                if (action_succeeded_) {
                    LiteVarAssignment lite_asgn;
                    FullAssignmentToLite(relax_asgn_, lite_asgn);
                    double asgn_rd = -1;
                    double asgn_rank = -1;
                    if (validator_.CheckRelaxation(lite_asgn, asgn_rd,
                                                   asgn_rank, true)) {
                        lite_asgn.clear();
                        FullAssignmentToLite(track_vars_asgn_, lite_asgn);
                        if (asgn_rd != -1) {
                            // multiply by 1000, since we don't accept doubles
                            lite_asgn.add_val(asgn_rd * 1000, 0);
                            if (asgn_rank != -1) {
                                lite_asgn.add_val(asgn_rank * 1000, 0);
                            }
                        }
                        validator_.sl_task_.ReportSolution(
                                asgns_.back().var_asgn_.mins_,
                                lite_asgn.mins_);
                    } else {
                        validator_.stats_.relax_post_filtered_.fetch_add(
                                1, std::memory_order_relaxed);
                        LOG4CXX_TRACE(logger, "Validator post-filtered: " <<
                                      asgns_.back());
                    }
                }
            } else if (last_action_.type_ == Action::Type::REMOTE_CHECK) {
                validator_.stats_.remote_check_.fetch_add(
                        1, std::memory_order_relaxed);
                LiteVarAssignment lite_asgn;
                FullAssignmentToLite(relax_asgn_, lite_asgn);
                double asgn_rd = -1;
                double asgn_rank = -1;
                if (action_succeeded_) {
                    const bool lrd_check = validator_.CheckRelaxation(
                            lite_asgn, asgn_rd, asgn_rank, true);
                    if (!lrd_check) {
                        validator_.stats_.relax_post_filtered_.fetch_add(
                                1, std::memory_order_relaxed);
                        action_succeeded_ = false;
                        LOG4CXX_TRACE(logger, "Validator post-filtered: " <<
                                      asgns_.back());
                    }
                    lite_asgn.clear();
                	FullAssignmentToLite(track_vars_asgn_, lite_asgn);
                	if (asgn_rd != -1) {
                        // multiply by 1000, since we don't accept doubles
                        lite_asgn.add_val(asgn_rd * 1000, 0);
                        if (asgn_rank != -1) {
                            lite_asgn.add_val(asgn_rank * 1000, 0);
                        }
                	}
                }
                validator_.SendForwardResult(last_action_.forward_id_,
                        action_succeeded_, lite_asgn.mins_);
            } else {
                // Simulation finished: we need to decide if we want to forward
                assert(action_succeeded_);
                adapter_->StopCollectingStats();
                CandidateAssignment &last_cas = asgns_.back();
                const bool forwarded =
                        validator_.CheckForward(
                                adapter_->GetCurrentStats().chunks_pos_,
                                last_cas);
                if (!forwarded) {
                    if (validator_.LocalZonesNumber() > 1) {
                        validator_.PushToLocalZone(
                                adapter_->GetCurrentStats().chunks_pos_,
                                last_cas);
                    } else {
                        // one zone, check here
                        last_cas.forw_id_ = -2; // switch to local check
                        asgns_.emplace_back(); // dummy (see pop_back() below)
                    }
                }
            }
            // Pop the assignment here, so we can use the info from it earlier
            asgns_.pop_back();
        }

        // check for searchlight termination
        if (validator_.CheckTerminate()) {
            // this will stop the search by failing the right branch
            LOG4CXX_ERROR(logger, "Terminating the validator by force");
            return solver->MakeFailDecision();
        }

        // pick the next assignment
        while (true) {
            if (asgns_.empty()) {
                if (!slave_) {
                    // need to ask the validator
                    CandidateVector *new_asgns =
                            validator_.GetNextAssignments();
                    if (!new_asgns) {
                        // No more assignments: stop the validator search
                        LOG4CXX_INFO(logger, "Stopping the validator search");
                        return solver->MakeFailDecision();
                    }
                    asgns_.swap(*new_asgns);
                    LOG4CXX_TRACE(logger, "Got " << asgns_.size()
                            << " new assignments to check");
                    delete new_asgns;
                } else {
                    LOG4CXX_INFO(logger, "Helper validator finished workload");
                    return solver->MakeFailDecision();
                }
            }
            assert(!asgns_.empty());
            const CandidateAssignment &ca = asgns_.back();
            if (validator_.CheckCandidateRD(ca)) {
                // Exit if this is not a relaxation or if it passes LRD
                break;
            } else {
                /*
                 * Didn't pass the LRD check -- ignore. We do the check based
                 * on the best RD obtained in the solver's leaf. Since we
                 * haven't had any new info since, the check is up-to-date.
                 */
                validator_.stats_.relax_pre_filtered_.fetch_add(
                        1, std::memory_order_relaxed);
                if (ca.forw_id_ >= 0) {
                    validator_.SendForwardResult(ca.forw_id_, false, {});
                }
                LOG4CXX_TRACE(logger, "Validator pre-filtered: " << ca);
                asgns_.pop_back();
            }
        }

        // Prepare assignment
        CandidateAssignment &ca = asgns_.back();
        /*
         * If we're relaxing, get the new maximum relaxation VC.
         * We need this, since candidates might fail the currently
         * non-relaxed constraint, but still qualify for a relaxed
         * result (by passing the LRD check).
         *
         * Note, this might modify the candidate's VC.
         */
        const int64_t forw_id = ca.forw_id_;
        LiteToFullAssignment(*last_asgn_, ca.var_asgn_);

        last_action_.forward_id_ = forw_id;
        if (forw_id >= 0) {
            last_action_.type_ = Action::Type::REMOTE_CHECK;
        } else if (forw_id == -2) {
            last_action_.type_ = Action::Type::LOCAL_CHECK;
        } else {
            last_action_.type_ =
                    validator_.forw_type_ != Validator::Forwarding::NONE ?
                    Action::Type::SIM : Action::Type::LOCAL_CHECK;
        }
        /*
         * If we're relaxing, compute maximum VC only when not simulating.
         */
        if (last_action_.type_ != Action::Type::SIM) {
            validator_.GetMaximumRelaxationVC(ca.relaxed_constrs_);
        }

        /*
         * Set adapter to either perform simulation or check immediately.
         */
        if (last_action_.type_ == Action::Type::SIM) {
            adapter_->SetAdapterMode(Adapter::DUMB);
            adapter_->StartCollectingStats();
        } else {
            adapter_->SetAdapterMode(Adapter::EXACT);
        }
        LOG4CXX_TRACE(logger, "Validator validating: " << ca);
        return solver->RevAlloc(new RestoreAssignment{last_asgn_.get(),
            ca.relaxed_constrs_, constrs_});
    }

    /**
     * Returns the name for debugging logging.
     * @return
     */
    virtual std::string DebugString() const {
        return std::string{"RestoreAssignment "} +
                (slave_ ? "(slave)" : "(master)");
    }

 private:
    /*
     * We need this Monitor to perform auxiliary tasks for the decision builder.
     * First of all, we need to restart the search to perform garbage
     * collecting. Secondly, we need to detect when the current assignment wa
     * validated successfully to produce NULL and, thus, signal a leaf to the
     * solver. This is done via AfterDecision().
     */
    class AuxRestoreMonitor : public SearchMonitor {
    public:
        /*
         * accept_flag is the flag to toggle at the decision builder after a
         * successful validation.
         */
        AuxRestoreMonitor(RestoreAssignmentBuilder &restore_db, Solver *solver,
                int restart_period) :
            SearchMonitor(solver),
            restore_db_(restore_db),
            restart_period_(restart_period) {}

        /*
         * Called after a decision was successfully accepted/refuted.
         */
        virtual void AfterDecision(Decision* const d, bool apply) {
            if (apply) {
                // Action was successful
                restore_db_.action_succeeded_ = true;
                if (restore_db_.last_action_.type_ != Action::Type::SIM) {
                	restore_db_.track_vars_asgn_.Store();
                	restore_db_.relax_asgn_.Store();
                }
                my_fail_ = true;
                solver()->Fail();
            } else {
                if (!my_fail_) {
                    // Action wasn't successful -- fail
                    restore_db_.action_succeeded_ = false;
                }
                // Restart the search for garbage collecting
                my_fail_ = false;
                if (restart_period_ &&
                        solver()->SearchDepth() > restart_period_) {
                    LOG4CXX_TRACE(logger,
                            "Restarting the validator for garbage collecting");
                    RestartCurrentSearch();
                    solver()->Fail();
                }
            }
        }

    private:
        // Flag to change when the next assignment is validated
        RestoreAssignmentBuilder &restore_db_;

        // Restart period of the validator
        int restart_period_;

        // True if we failed the accept branch
        bool my_fail_ = false;
    };

    // The servicing validator. For slaves -- the parent validator.
    Validator &validator_;

    // Adapter used for access
    AdapterPtr adapter_;

    // Relaxable constraints
    const RelaxableConstraints &constrs_;

    // Assignments to validate
    CandidateVector asgns_;

    // True, if the last action was a success (no fail on solver)
    bool action_succeeded_;

    // The auxiliary monitor
    AuxRestoreMonitor aux_monitor_;

    // Assignment used in the last action
    AssignmentPtr last_asgn_;

    // Additional tracking vars
    Assignment track_vars_asgn_;

    // Relaxable constraints expressions prototype
    Assignment relax_asgn_;

    // Last action taken by the DB
    Action last_action_;

    // True if this is a slave builder
    const bool slave_;
};

Validator::Validator(Searchlight &sl, SearchlightTask &sl_task,
        const StringVector &var_names) :
        sl_(sl),
        sl_task_(sl_task),
        relaxator_(sl.GetRelaxator()),
        solver_("validator solver"),
        adapter_(sl.CreateAdapter("validator")), // DUMB mode by default!
        search_vars_prototype_(&solver_),
		track_vars_prototype_(&solver_),
		rel_const_prototype_(&solver_) {

    // First, clone the solver (solver 0 is always available)
    ExportModel(sl_.GetSearchSolver(0), &validator_model_);
    if (!CloneModel(sl_, validator_model_, solver_, adapter_)) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Cannot create the validator solver!";
    }

    // logging
    if (logger->isDebugEnabled()) {
        /*
         *  Unfortunately, it will dump it into std::cerr, which is
         *  redirected by SciDb to a file.
         */
        std::cerr << "Logging model for validator...\n";
        ModelVisitor *pmv = solver_.MakePrintModelVisitor();
        solver_.Accept(pmv);
    }

    /*
     * Now, we want to find all decision variables. We do this via a
     * custom visitor.
     */
    ModelCollector var_finder;
    solver_.Accept(&var_finder);
    var_finder.FillAssignment(var_names, search_vars_prototype_);
    var_finder.FillAssignment(sl.GetTrackExprs(), track_vars_prototype_);
    relaxable_constrs_ = var_finder.GetRelaxableConstraints();
    var_finder.FillRelaxAssignment(rel_const_prototype_);

    const SearchlightConfig &sl_config = sl_task_.GetConfig();
    max_pending_validations_ =
            sl_config.get("searchlight.validator.max_validations", 1000000);
    restart_period_ = sl_config.get("searchlight.validator.restart_period",
            1024);
    max_helpers_allowed_ =
            sl_config.get("searchlight.validator.max_helpers", 1);
    helper_workload_ =
            sl_config.get("searchlight.validator.helper_workload", 1000);
    send_info_period_ =
            sl_config.get("searchlight.validator.send_info_period", 0);

    // Are we using dynamic scheduling for helpers?
    dynamic_helper_scheduling_ =
            sl_config.get("searchlight.dynamic_scheduling", 1);

    // Watermarks
    low_watermark_ = sl_config.get("searchlight.validator.low_watermark", 500);
    high_watermark_ =
            sl_config.get("searchlight.validator.high_watermark", 1000);
    rebal_watermark_ =
            sl_config.get("searchlight.validator.rebal_watermark", 10000);

    // Determine forwarding type
    const std::string forw_type =
            sl_config.get("balance.validator_balance", "none");
    if (forw_type == "stripes") {
        forw_type_ = Forwarding::STRIPES;
        LOG4CXX_INFO(logger, "Forwarding is stripe-based...");
    } else if (forw_type == "dynamic") {
        forw_type_ = Forwarding::DYNAMIC;
        LOG4CXX_INFO(logger, "Forwarding is dynamic mapping...");
    } else {
        forw_type_ = Forwarding::NONE;
        LOG4CXX_INFO(logger, "Forwarding is disabled for the validator...");
    }

    // Determine logical id
    const auto &active_validators = sl_task.GetActiveValidators();
    const auto iter = std::find(active_validators.begin(),
            active_validators.end(), sl_task.GetInstanceID());
    my_logical_id_ = iter == active_validators.end() ? -1 :
            iter - active_validators.begin();

    // Configured number of zones
    size_t zones_num = sl_config.get("searchlight.validator.zones", 0);
    const auto &search_array_desc = adapter_->GetSearchArrayDesc();

    // If not specified, determine automatically
    if (zones_num == 0) {
        /*
         * Determine the number of zones. We want at least two zones to fit in
         * memory. However, we should take into account two types of memory:
         * local (L) for disk buffering and temporary (T) for temporary
         * transferred arrays.
         */
        const unsigned int ZONES_IN_MEMORY = 2;

        // Params
        const size_t validator_size =
                (search_array_desc.GetSearchArraySize() - 1) /
                active_validators.size() + 1;
        const size_t cache_size = size_t(Config::getInstance()->getOption<int>(
                scidb::CONFIG_SMGR_CACHE_SIZE)) * 1024 * 1024;
        const size_t tmp_size = size_t(Config::getInstance()->getOption<size_t>(
                scidb::CONFIG_MEM_ARRAY_THRESHOLD)) * 1024 * 1024;
        const size_t inst_count = sl_task_.GetQueryInstanceCount();

        // Number of zones
        zones_num = std::max(ZONES_IN_MEMORY * validator_size /
                inst_count / cache_size,
                ZONES_IN_MEMORY * validator_size * (inst_count - 1) /
                inst_count / tmp_size) + 1;
    }

    /*
     * Candidate sorting order.
     */
    ValidatorQueue::Ordering ca_ordering = ValidatorQueue::Ordering::COST;
    const std::string cand_sort = sl_config.get("searchlight.validator.sort",
            "cost");
    if (cand_sort == "rd") {
        ca_ordering = ValidatorQueue::Ordering::RD;
        // FIXME: temporary turn off zones for RD
        zones_num = 1;
        LOG4CXX_WARN(logger, "RD sorting forces a single validator zone!");
    } else if (cand_sort != "cost") {
        LOG4CXX_WARN(logger, "Unknown candidate ordering. Defaulting to cost.");
    }

    /*
     *  Init structures.
     *
     *  FIXME: We do it even for inactive validators despite the fact they
     *  don't have zones. We do it to make some functions work, like
     *  LocalZonesNumber(). Shoul revisit this.
     */
    val_queue_.reset(new ValidatorQueue(zones_num, ca_ordering, my_logical_id_,
            active_validators.size()));
    LOG4CXX_INFO(logger, "Creating " << zones_num
            << " zones for the validator...");
    // Need to know the zones for each validator
    for (size_t val_id = 0; val_id < active_validators.size(); ++val_id) {
        chunk_zones_.emplace_back(search_array_desc.CreateChunkZones(
                {active_validators.size(), zones_num},
                {val_id}));
    }
}

Validator::~Validator() {
    std::ostringstream os;
    os << stats_;
    LOG4CXX_INFO(logger, os.str());
}

void Validator::AddSolution(CandidateAssignment &&cas) {
    std::unique_lock<std::mutex> validate_lock(to_validate_mtx_);
    stats_.total_local_cands_.fetch_add(1, std::memory_order_relaxed);
    cas.forw_id_ = -1;
    LOG4CXX_TRACE(logger, "Added local candidate: " << cas);
    PushCandidate(std::move(cas), 0 /* doesn't matter */);

    // We might have to wait -- flood control
    while (val_queue_->TotalCands() > max_pending_validations_) {
        /*
         * It is safe to use the same condition since the validator will be
         * certainly non-blocked.
         */
        LOG4CXX_DEBUG(logger, "Waiting for the validator to catch up"
                ", queue size=" << val_queue_->TotalCands());
        validate_cond_.wait(validate_lock);
    }

    // Notify the validator in case it is blocked and unlock
    validate_cond_.notify_all();
}

void Validator::PushCandidate(CandidateAssignment &&asgn, size_t zone) {
    const bool non_sim_cand = asgn.forw_id_ == -1;
    val_queue_->PushCandidate(std::move(asgn), zone);
    if (my_logical_id_ != -1 && !non_sim_cand) {
        BroadcastCandidatesCount();
    }
}

void Validator::BroadcastCandidatesCount() {
    if (my_logical_id_ != -1 && send_info_period_ > 0) {
        const size_t local_cands = val_queue_->ReadyCandsNum();
        // Active validator: check if we want to send the info update
        const size_t cands_diff = std::abs(int64_t(local_cands) -
                int64_t(last_sent_cands_num_));
        if (cands_diff >= send_info_period_ ||
                (local_cands == 0 && last_sent_cands_num_ != 0)) {
            // Send if the diff is high or we're finished
            last_sent_cands_num_ = local_cands;
            sl_task_.BroadcastValidatorInfo(local_cands);
        }
    }
}

void Validator::AddRemoteCandidates(CandidateVector &cands,
        const Int64Vector &zones,
        InstanceID src, uint64_t forw_id) {
    assert(cands.size() == zones.size());

    std::lock_guard<std::mutex> validate_lock(to_validate_mtx_);
    for (size_t i = 0; i < cands.size(); i++) {
        const int64_t cand_id = remote_cand_id_++;
        cands[i].forw_id_ = cand_id;
        LOG4CXX_TRACE(logger, "Added remote candidate: " << cands[i]);
        PushCandidate(std::move(cands[i]), zones[i]);
        remote_candidates_.emplace(cand_id, std::make_pair(src, forw_id++));
    }
    stats_.total_remote_cands_.fetch_add(
            cands.size(), std::memory_order_relaxed);
    validate_cond_.notify_all();
}

void Validator::HandleForwardResult(uint64_t id, bool result,
		const std::vector<int64_t> &add_vals) {
    if (IsReforward(id)) {
        // Forward to the original sender
        std::unique_lock<std::mutex> validate_lock{to_validate_mtx_};
        assert(reforwarded_candidates_.find(id) !=
                reforwarded_candidates_.end());
        const int64_t original_id = reforwarded_candidates_[id];
        reforwarded_candidates_.erase(id);
        if (reforwarded_candidates_.empty()) {
            validate_cond_.notify_all();
        }
        validate_lock.unlock();
        SendForwardResult(original_id, result, add_vals);
    } else {
        std::lock_guard<std::mutex> validate_lock{to_validate_mtx_};
        assert(forwarded_candidates_.find(id) != forwarded_candidates_.end());
        if (result) {
            sl_task_.ReportSolution(forwarded_candidates_[id].mins_, add_vals);
        }
        forwarded_candidates_.erase(id);
        if (forwarded_candidates_.empty()) {
            // Check if we went idle
            if (Idle_NL()) {
                sl_.EventNotify(SearchlightEvent::VALIDATOR_IDLE);
            }
            validate_cond_.notify_all();
        }
    }
}

void Validator::SendForwardResult(int64_t forw_id, bool result,
		const std::vector<int64_t> &add_vals) {
    assert(forw_id >= 0);
    std::unique_lock<std::mutex> validate_lock(to_validate_mtx_);
    const auto forw_info = remote_candidates_[forw_id];
    remote_candidates_.erase(forw_id);
    validate_lock.unlock();

    sl_task_.SendBalanceResult(forw_info.first, forw_info.second, result,
    		add_vals);
}

template <typename CoordinatesSequence>
int Validator::DetermineLocalZone(const CoordinatesSequence &chunks,
        size_t val) const {
    if (LocalZonesNumber() == 1) {
        return 0;
    }
    std::vector<int> zone_counts(LocalZonesNumber()); // zones number
    size_t max_zone = 0;
    // zones[1] contains local zones; zones[0] -- global ones
    adapter_->GetSearchArrayDesc().GetStripesChunkDistribution(
            chunks, zone_counts, max_zone, chunk_zones_[val].zones_[1]);

    // Narrowing cast here, since zone number are ints.
    return static_cast<int>(max_zone);
}

void Validator::PushToLocalZone(const CoordinateSet &chunks,
                                const CandidateAssignment &asgn) {
    // Determine zone
    assert(my_logical_id_ != -1);
    const size_t zone = DetermineLocalZone(chunks, my_logical_id_);
    CandidateAssignment new_asgn(asgn);
    new_asgn.forw_id_ = -2;
    std::lock_guard<std::mutex> lock{to_validate_mtx_};
    PushCandidate(std::move(new_asgn), zone);
}

bool Validator::CheckForward(const CoordinateSet &chunks,
                             const CandidateAssignment &asgn) {
    const auto &active_validators = sl_task_.GetActiveValidators();
    std::vector<int> inst_counts(active_validators.size());

    // Get the distribution
    size_t max_inst = 0;
    switch (forw_type_) {
        case Forwarding::STRIPES:
            // zones_[0] contains global validator zones
            // doesn't matter which spec to use: global ones are the same
            adapter_->GetSearchArrayDesc().GetStripesChunkDistribution(
                    chunks, inst_counts, max_inst, chunk_zones_[0].zones_[0]);
            break;
        case Forwarding::DYNAMIC: {
            std::vector<int> inst_counts_all(sl_task_.GetQueryInstanceCount());
            adapter_->GetSearchArrayDesc().GetDynamicChunksDistribution(
                    sl_task_.GetQueryContext(), chunks, inst_counts_all);
            for (size_t i = 0; i < active_validators.size(); i++) {
                inst_counts[i] = inst_counts_all[active_validators[i]];
            }
            max_inst = std::distance(inst_counts.begin(),
                    std::max_element(inst_counts.begin(), inst_counts.end()));
            break;
        }
        case Forwarding::NONE:
            assert(false);
            break;
    }

    // A couple of forward heuristsics...
    if (my_logical_id_ != -1 && max_inst != my_logical_id_ &&
            inst_counts[my_logical_id_] == inst_counts[max_inst]) {
        /*
         *  In case the local instance contains the same number of chunks,
         *  we don't forward.
         */
        max_inst = my_logical_id_;
    }

    /*
     * Check if we have a borderline tie. In which case we use the
     * validator's info.
     */
    int tie_inst = 0;
    if (max_inst > 0 && inst_counts[max_inst - 1] ==
            inst_counts[max_inst]) {
        tie_inst = -1;
    } else if (max_inst < inst_counts.size() - 1 &&
            inst_counts[max_inst + 1] == inst_counts[max_inst]) {
        tie_inst = 1;
    }

    const auto &val_cands_info = val_queue_->ValidatorsCandsInfo();
    // FIXME: dangerous code, since val_cands_info is non-atomic!!!
    if (val_cands_info[max_inst + tie_inst] < val_cands_info[max_inst]) {
        // Tie-break in favor of the instance with less candidates
        LOG4CXX_TRACE(logger, "Tie-break of inst=" << max_inst <<
                " in favor of inst=" << max_inst + tie_inst);
        max_inst += tie_inst;
    }

    if (my_logical_id_ == -1 || max_inst != my_logical_id_) {
        // Update stats locally (global will follow through broadcasts)
        UpdateCandsInfo(val_cands_info[max_inst] + 1, max_inst);

        // forward
        uint64_t cand_id;
        {
            std::lock_guard<std::mutex> lock{to_validate_mtx_};
            cand_id = forw_id_++;
            forwarded_candidates_.emplace(cand_id, asgn.var_asgn_);
        }

        /*
         * Compute the candidate's zone. Use int64 to conform
         * to the message's format.
         */
        const Int64Vector zones{DetermineLocalZone(chunks, max_inst)};

        // Prepare coordinates
        // FIXME: I'm leaving this code for now for future reference (AK)
//        std::vector<std::vector<int64_t>> coords;
//        if (LocalZonesNumber() > 1) {
//            // Need coordinate de-duplication
//            // Don't send coordinates in case of the single zone
//            std::unordered_set<int64_t> dedup_coords;
//            const size_t coord_ord = chunk_zones_.zones_[1].dim_num_;
//            for (const auto &c: chunks) {
//                dedup_coords.insert(c[coord_ord]);
//            }
//            coords.emplace_back(dedup_coords.begin(), dedup_coords.end());
//        }
        stats_.total_forw_.fetch_add(1, std::memory_order_relaxed);
        sl_task_.ForwardCandidates({asgn}, zones, active_validators[max_inst],
                cand_id);
        return true;
    }

    // This validator must be active here!
    assert(my_logical_id_ != -1);
    return false;
}

bool Validator::CheckRelaxation(const LiteVarAssignment &relax_asgn, double &rd,
                                double &rank, bool report_rd) const {
    if (!relaxator_) {
        // Not relaxing; check always passes
        return true;
    }
    return relaxator_->CheckResult(relax_asgn, report_rd, rd, rank);
}

bool Validator::CheckCandidateRD(const CandidateAssignment &ca) const {
    assert(ca.relaxed_constrs_.empty() || relaxator_);
    if (!relaxator_) return true;
    // First check the relaxing
    bool res = ca.relaxed_constrs_.empty() ||
            ca.best_rd_ <= relaxator_->GetLRD();
    // Check the contractor if we're not relaxing anymore
    if (!relaxator_->Relaxing()) {
        if (!ca.rc_vals_.empty()) {
            res = relaxator_->CheckContraction(ca.rc_vals_);
        } else {
            res = relaxator_->CheckContractionRank(ca.best_rank_);
        }
    }
    return res;
}

void Validator::GetMaximumRelaxationVC(Int64Vector &vc) const {
    if (relaxator_) {
        assert(vc.size() % 3 == 0);
        // vc spec consists of triplets
        const size_t viol_num = vc.size() / 3;
        /*
         * Note, viol_num might be 0, if there're no violations. In this case,
         * the validated result might still be a relaxed one and pass LRD.
         * We use 1 probable violation to determine the maximum intervals.
         */
        vc = relaxator_->GetMaximumRelaxationVCSpec(viol_num ? viol_num : 1);
    }
}

int Validator::FindValidatorForReforwards() {
    assert(my_logical_id_ >= 0);

    // Threshold for determining an idle validator (some small value)
    const size_t BUSY_THRESHOLD = 100;

    // First, check recent validators
    // FIXME: dangerous code, since val_cands_info is non-atomic!!!
    const auto &val_cands_info = val_queue_->ValidatorsCandsInfo();
    for (auto rit = validators_mru_reforw_.rbegin();
            rit != validators_mru_reforw_.rend(); ++rit) {
        const int val_id = *rit;
        if (val_cands_info[val_id] < BUSY_THRESHOLD) {
            if (rit != validators_mru_reforw_.rbegin()) {
                validators_mru_reforw_.erase((++rit).base());
                validators_mru_reforw_.push_back(val_id);
            }
            return val_id;
        }
    }

    // Then, check left and right neighbors (since they overlap)
    if (my_logical_id_ > 0) {
        // left
        const int val_id = my_logical_id_ - 1;
        if (val_cands_info[val_id] < BUSY_THRESHOLD) {
            validators_mru_reforw_.push_back(val_id);
            return val_id;
        }
    }
    if (my_logical_id_ < val_cands_info.size() - 1) {
        // right
        const int val_id = my_logical_id_ + 1;
        if (val_cands_info[val_id] < BUSY_THRESHOLD) {
            validators_mru_reforw_.push_back(val_id);
            return val_id;
        }
    }

    // Then, just try to find anybody...
    for (int i = 0; i < val_cands_info.size(); ++i) {
        if (i != my_logical_id_ && val_cands_info[i] < BUSY_THRESHOLD) {
            validators_mru_reforw_.push_back(i);
            return i;
        }
    }

    // No idle validators
    return -1;
}

void Validator::Synchronize() const {
    std::unique_lock<std::mutex> validate_lock(to_validate_mtx_);
    while (!val_queue_->Empty()) {
        LOG4CXX_DEBUG(logger, "Synchronizing with the validator"
                ", queue size=" << val_queue_->TotalCands());
        validate_cond_.wait(validate_lock);
    }
}

void Validator::operator()() {
    LOG4CXX_INFO(logger, "Starting the validator search");
    DecisionBuilder *db =
            solver_.RevAlloc(new RestoreAssignmentBuilder(*this, adapter_,
                    &solver_, relaxable_constrs_,
                    restart_period_, &search_vars_prototype_,
					track_vars_prototype_, rel_const_prototype_, {},
                    false));
    try {
        solver_status_ = solver_.Solve(db);
    } catch (const scidb::Exception &ex) {
        // Might catch an exception from SearchlightTask
        LOG4CXX_ERROR(logger, "Validator caught an exception");
        sl_task_.HandleSearchlightError(ex.copy());
        // Handle possible waiters
        GetNextAssignments();
    }

    // Check if Solve() ended because of the false model
    const Searchlight::Status sl_status = sl_.GetStatus();
    if (sl_status != Searchlight::Status::TERMINATED
            && sl_status != Searchlight::Status::COMMITTED) {
        /*
         * Solve() exited because of the model. We still have to handle the
         * main loop and all control messages.
         */
        GetNextAssignments();
    }
    // Should wait for helpers
    if (!validator_helpers_.empty()) {
        LOG4CXX_INFO(logger, "Validator waiting for helpers to finish...");
        validator_helpers_.clear();
    }
    LOG4CXX_INFO(logger, "Validator exited solve()");
}

CandidateVector *Validator::GetNextAssignments() {
    /*
     * For now, a simple policy: wait until we get an assignment to
     * validate and then check it.
     */
    std::unique_lock<std::mutex> validate_lock(to_validate_mtx_);
    while (true) {
        // Check interesting statuses first
        const Searchlight::Status sl_status = sl_.GetStatus();
        if (sl_status == Searchlight::Status::TERMINATED) {
            LOG4CXX_ERROR(logger, "Terminating validator by force");
            // Notify the main solver about catching up, if needed
            val_queue_->Clear();
            validate_cond_.notify_all();
            return nullptr;
        } else if (sl_status == Searchlight::Status::COMMITTED) {
            assert(FinishedLocally());
            LOG4CXX_INFO(logger, "Committing validator");
            return nullptr;
        } else if (sl_status == Searchlight::Status::FIN_SEARCH) {
            if (FinishedLocally()) {
                /*
                 * Have to avoid recursive locking in case we're at the
                 * coordinator and it immediately "sends" commit.
                 *
                 * Since after re-locking the status might have changed,
                 * re-check it again.
                 *
                 * Note, we still cannot exit the loop: might have forwards
                 * coming up.
                 */
                validate_lock.unlock();
                LOG4CXX_INFO(logger, "Validator finished locally");
                sl_.ReportFinValidator();
                validate_lock.lock();
                continue;
            }
        }
        // Then check if we have any candidates or forwards
        if (!val_queue_->Empty()) {
            break;
        }
        // Nothing new to do -- wait
        validate_cond_.wait(validate_lock);
    }

    // This is mine!
    assert(!val_queue_->Empty());
    CandidateVector *res = new CandidateVector;
    *res = GetWorkload();

    // Check if we need to rebalance
    const size_t cands_ready = val_queue_->ReadyCandsNum();
    if (cands_ready > rebal_watermark_) {
        const int val_rebal = FindValidatorForReforwards();
        if (val_rebal >= 0) {
            /*
             * How many candidates? Rule-of-thumb. We don't want to overload
             * the other validator and at the same time leave this one with
             * too few candidates.
             */
            const size_t cands_rebal = std::min(
                    cands_ready / 2,
                    rebal_watermark_);
            RebalanceAndSend(cands_rebal, val_rebal);
        }
    }

    // Try to off-load the rest to helpers
    bool should_run_helpers = !val_queue_->Empty();
    while (should_run_helpers) {
        should_run_helpers = false;
        if (sl_task_.ReserveThread()) {
           bool pers_helper;
           should_run_helpers = AddValidatorHelperInt(&pers_helper) &&
                   !val_queue_->Empty();
        }
    }

    // Check if we went idle
    if (Idle_NL()) {
        sl_.EventNotify(SearchlightEvent::VALIDATOR_IDLE);
    }

    // Flow control: notify the main solver about catching up
    validate_cond_.notify_all();

    return res;
}

void Validator::RebalanceAndSend(size_t cands, int validator) {
    // Sanity checks
    assert(cands > 0);
    const auto &active_validators = sl_task_.GetActiveValidators();
    assert(validator >= 0 && validator < active_validators.size() &&
            validator != my_logical_id_);

    /*
     * Prepare info. The main reason we need two lists is that we have
     * different id ranges for forwards and re-forwards. So, we send them
     * as two different batches.
     */
    CandidateVector reforwards, forwards;
    Int64Vector reforwards_zones, forwards_zones;
    reforwards.reserve(cands);
    forwards.reserve(cands);
    reforwards_zones.reserve(cands);
    forwards_zones.reserve(cands);
    const uint64_t start_forw_id = forw_id_;
    const uint64_t start_reforw_id = reforw_id_;

    // Rebalance and fill in forward/reforward info
    val_queue_->FillInRebalanceTransfer(cands, reforwards, reforwards_zones,
            forwards, forwards_zones);
    assert(!reforwards.empty() || !forwards.empty());
    for (const auto &ca: reforwards) {
        reforwarded_candidates_.emplace(reforw_id_++, ca.forw_id_);
    }
    for (const auto &ca: forwards) {
        forwarded_candidates_.emplace(forw_id_++, ca.var_asgn_);
    }

    // Send
    LOG4CXX_INFO(logger, "Rebalancing: sending " <<
            reforwards.size() + forwards.size() <<
            " candidates to validator " << validator);
    if (!reforwards.empty()) {
        sl_task_.ForwardCandidates(reforwards, reforwards_zones,
                active_validators[validator], start_reforw_id);
        stats_.total_reforw_.fetch_add(
                reforwards.size(), std::memory_order_relaxed);
    }
    if (!forwards.empty()) {
        sl_task_.ForwardCandidates(forwards, forwards_zones,
                active_validators[validator], start_forw_id);
        stats_.total_forw_.fetch_add(
                forwards.size(), std::memory_order_relaxed);
    }
}

void ValidatorQueue::FillInRebalanceTransfer(size_t cands,
                                             CandidateVector &reforwards,
                                             Int64Vector &reforwards_zones,
                                             CandidateVector &forwards,
                                             Int64Vector &forwards_zones) {
    // Sanity checks
    assert(cands < ReadyCandsNum() && val_logical_id_ != -1);
    // We go through zones in the LRU order (from old to new ones)
    for (auto zone: zones_mru_) {
        auto &zone_cands = zones_[zone];
        while (cands && !zone_cands.empty()) {
            CandidateAssignment &ca = const_cast<CandidateAssignment &>(
                    zone_cands.top());
            if (ca.forw_id_ >= 0) {
                // Reforward
                reforwards.push_back(std::move(ca));
                reforwards_zones.push_back(zone);
            } else {
                // Local
                assert(ca.forw_id_ == -2);
                forwards.push_back(std::move(ca));
                forwards_zones.push_back(zone);
            }
            zone_cands.pop();
            --cands;
        }
        if (!cands) break;
    }

    // Accounting
    const size_t total_rebalanced = reforwards.size() + forwards.size();
    assert(total_rebalanced < to_validate_total_);
    to_validate_total_.fetch_sub(total_rebalanced, std::memory_order_relaxed);
    assert(total_rebalanced < validators_cands_info_[val_logical_id_]);
    validators_cands_info_[val_logical_id_] -= total_rebalanced;
}

CandidateVector ValidatorQueue::GetCandidates(size_t num) {
    assert(to_validate_total_.load(std::memory_order_relaxed));
    CandidateVector res;
    res.reserve(num);
    // First try to find some ready candidates
    if (ReadyCandsNum() > 0) {
        // Determine the zone
        size_t zone = 0;
        if (zones_mru_.size() > 1) {
            /*
             * Simple cost-heuristic. Possibly switch to a zone with better
             * candidates, but still keep in the same zone.
             *
             * TODO: revisit later. Something more sophisticated, like taking
             * candidates from different zones and then batching them together
             * to avoid thrashing between zones.
             */
            auto rit = zones_mru_.rbegin();
            // Find non-empty zone: note, should have at least one
            while (zones_[*rit].empty()) {
                ++rit;
                assert(rit != zones_mru_.rend());
            }
            auto zone_rit = rit;
            double min_cost = CandidateCost(zones_[*rit].top());
            while (++rit != zones_mru_.rend()) {
                const double new_cost = CandidateCost(zones_[*rit].top());
                if (new_cost < min_cost &&
                        (min_cost - new_cost) / min_cost > 0.1) {
                    // Some "significant" change
                    min_cost = new_cost;
                    zone_rit = rit;
                }
            }
            // zone_rit points to the minimum
            zone = *zone_rit;
            if (zone_rit != zones_mru_.rbegin()) {
                // Put the zone as first in the MRU
                zones_mru_.erase((++zone_rit).base());
                zones_mru_.push_back(zone);
            }
        }
        // We have the zone now
        auto &zone_queue = zones_[zone];
        while (res.size() < num && !zone_queue.empty()) {
            // const_cast is okay here for move
            res.push_back(std::move(
                    const_cast<CandidateAssignment &>(zone_queue.top())));
            zone_queue.pop();
        }
        // Update cands info
        assert(val_logical_id_ != -1);
        validators_cands_info_[val_logical_id_] -= res.size();
    }
    /*
     *  Now about non-simulated candidates. Take 'num' of them as well.
     */
    while (num > 0 && !nonsim_cands_.empty()) {
        res.push_back(std::move(nonsim_cands_.front()));
        nonsim_cands_.pop();
        --num;
    }
    // Accounting
    assert(!res.empty());
    to_validate_total_.fetch_sub(res.size(), std::memory_order_relaxed);
    return res;
}

CandidateVector Validator::GetWorkload() {
    CandidateVector res = val_queue_->GetCandidates(helper_workload_);
    if (logger->isDebugEnabled() && !val_queue_->Empty()) {
        LOG4CXX_DEBUG(logger, "Dispatched a new workload, left=" <<
            val_queue_->TotalCands());
    }
    BroadcastCandidatesCount();
    return res;
}

bool Validator::AddValidatorHelper(bool *persistent) {
    std::lock_guard<std::mutex> validate_lock{to_validate_mtx_};
    return AddValidatorHelperInt(persistent);
}

bool Validator::AddValidatorHelperInt(bool *persistent) {
    Validator::ValidatorHelper *helper = nullptr;
    if (!val_queue_->Empty()) {
        helper = DispatchValidatorHelper();
    }
    if (!helper) {
        if (dynamic_helper_scheduling_) {
            LOG4CXX_DEBUG(logger, "No helper needed, returning the thread...");
            sl_task_.FreeThread();
        }
        return false;
    }

    LOG4CXX_INFO(logger, "Added another validator helper...");
    if (val_queue_->TotalCands() > high_watermark_) {
        helper->SetPersistent(true);
        *persistent = true;
        LOG4CXX_INFO(logger, "The helper is made persistent...");
    } else {
        *persistent = false;
    }

    helper->RunWorkload(GetWorkload());
    return true;
}

Validator::ValidatorHelper *Validator::DispatchValidatorHelper() {
    if (!free_validator_helpers_.empty()) {
        const int id = free_validator_helpers_.front();
        free_validator_helpers_.pop_front();
        ValidatorHelper *helper = validator_helpers_[id].get();
        return helper;
    } else if (dynamic_helper_scheduling_ ||
            validator_helpers_.size() < max_helpers_allowed_) {
        const int id = validator_helpers_.size();
        ValidatorHelper *helper = nullptr;
        try {
            helper = new ValidatorHelper{id, *this, search_vars_prototype_};
            validator_helpers_.emplace(id, ValidatorHelper::UniquePtr{helper});
        } catch (const scidb::Exception &) {
            // Cannot create a helper. This is not fatal.
            LOG4CXX_ERROR(logger, "Cannot create a validator helper. "
                    "Possibly an error during cloning the model.");
        } catch (...) {
            // Cannot create a helper. This is not fatal.
            LOG4CXX_ERROR(logger, "Cannot create a validator helper. "
                    "Reason unknown.");
        }
        return helper;
    }

    return nullptr;
}

Validator::ValidatorHelper::ValidatorHelper(int id, Validator &parent,
        const Assignment &init_assignment) :
                parent_(parent),
                solver_{"validator helper solver " + std::to_string(id)},
                adapter_{parent_.sl_.CreateAdapter(
                        "validator helper_" + std::to_string(id))},
                prototype_{&solver_},
				track_prototype_{&solver_},
				rel_const_prototype_{&solver_},
                id_{id} {
    /*
     *  First, clone the solver. It's crucial to create it from the model and
     *  not from the validator's solver. The validator solver might be
     *  running when the helper is being created, and, thus, the variables
     *  might have very different domains.
     */
    if (!CloneModel(parent_.sl_, parent_.validator_model_, solver_, adapter_)) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Cannot create a validator helper solver!";
    }

    /*
     * Now, we want to find all decision variables. We do this via a
     * custom visitor.
     */
    ModelCollector var_finder;
    solver_.Accept(&var_finder);
    StringVector search_var_names;
    for (const auto &elem: init_assignment.IntVarContainer().elements()) {
    	search_var_names.push_back(elem.Var()->name());
    }
    var_finder.FillAssignment(search_var_names, prototype_);
    var_finder.FillAssignment(parent_.sl_.GetTrackExprs(), track_prototype_);
    relaxable_constrs_ = var_finder.GetRelaxableConstraints();
    var_finder.FillRelaxAssignment(rel_const_prototype_);
}

void Validator::ValidatorHelper::operator()() {
    do {
        LOG4CXX_INFO(logger, "Starting a helper validator workload...");
        /*
         * Using uniue_ptr here since we create multiple dbs that are not
         * going to be freed by the solver -- they will reside before the
         * search sentinel.
         */
        const std::unique_ptr<DecisionBuilder> db{
                new RestoreAssignmentBuilder(parent_, adapter_,
                        &solver_, relaxable_constrs_,
                        0, &prototype_, track_prototype_,
                        rel_const_prototype_,
                        std::move(workload_), true)};
        adapter_->SetAdapterMode(Adapter::DUMB); // Will be switched in Next()
        try {
            solver_.Solve(db.get());
        } catch (const scidb::Exception &ex) {
            // Might catch an exception from SearchlightTask
            parent_.sl_task_.HandleSearchlightError(ex.copy());
            LOG4CXX_ERROR(logger, "Helper caught an exception");
            break;
        }
        workload_.clear();
        LOG4CXX_INFO(logger, "Validator helper finished workload");

        // Need to determine if we want to continue
        bool should_stop = false;
        if (parent_.CheckTerminate()) {
            LOG4CXX_ERROR(logger, "Terminating validator helper by force!");
            break;
        }
        if (parent_.dynamic_helper_scheduling_) {
            const size_t cands_num =
                    parent_.val_queue_->TotalCands();
            if (persistent_) {
                /*
                 * If we're persistent, we will notify the task only when we've
                 * fallen below the low watermark.
                 */
                if (cands_num < parent_.low_watermark_) {
                    persistent_ = false;
                    should_stop = parent_.sl_task_.PendingSolverJobs(true);
                    LOG4CXX_INFO(logger, "Persistent helper is below the "
                            "low watermark. Switching off persistence...");
                }
            } else {
                /*
                 * Non-persistent. We check every time we're below the high
                 * watermark if the solver needs threads.
                 */
                if (cands_num < parent_.high_watermark_) {
                    should_stop = parent_.sl_task_.PendingSolverJobs(false);
                }
            }
        }

        std::lock_guard<std::mutex> lock{parent_.to_validate_mtx_};
        if (!should_stop && !parent_.val_queue_->Empty()) {
            SetWorkload(parent_.GetWorkload());
        }
    } while (!workload_.empty());

    std::lock_guard<std::mutex> lock{parent_.to_validate_mtx_};
    parent_.free_validator_helpers_.push_back(id_);
    parent_.validate_cond_.notify_all();
    parent_.sl_task_.FreeThread();
    persistent_ = false;

    LOG4CXX_DEBUG(logger, "Helper finished, total assignments=" <<
            parent_.val_queue_->TotalCands());
    LOG4CXX_INFO(logger, "Validator helper exited solve()");
}

void Validator::ValidatorHelper::RunWorkload(CandidateVector &&workload) {
    assert(workload_.empty());
    Prepare();
    SetWorkload(std::move(workload));
    thr_ = std::thread{std::ref(*this)};
}

} /* namespace searchlight */

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
#include "searchlight_collector.h"
#include "searchlight_task.h"

#include "ortools_model.h"

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.validator"));

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
     */
    RestoreAssignment(Assignment *asgn) : asgn_(asgn) {}

    /**
     * The left branch. In this case it restores the assignment.
     *
     * @param s the solver
     */
    virtual void Apply(Solver* const s) {
        LOG4CXX_TRACE(logger, "Validating: " << asgn_->DebugString());
        asgn_->Restore();
    }

    /**
     * The right branch. It does nothing.
     *
     * @param s the solver
     */
    virtual void Refute(Solver* const s) {}

    /**
     * Returns a string for debug logging.
     *
     * @return a debug string
     */
    virtual std::string DebugString() const {
        return "RestoreAssignmentDecision";
    }

private:
    // The assignment to restore
    Assignment *asgn_;
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
        int forward_id_;
    };
public:

    /**
     * Creates a new restore assignment builder.
     *
     * @param validator the validator producing assignments
     * @param s the solver
     * @param restart_period validator restart period
     * @param prototype prototype assignment to generate validations
     */
    explicit RestoreAssignmentBuilder(Validator &validator, Solver *s,
            int restart_period, const Assignment *prototype)
        : validator_(validator),
          action_succeeded_(false),
          aux_monitor_(*this, s, restart_period),
          last_asgn_{new Assignment{prototype}} {}


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
                if (action_succeeded_) {
                    LiteVarAssignment lite_asgn;
                    FullAssignmentToLite(*last_asgn_, lite_asgn);
                    validator_.sl_task_.ReportSolution(lite_asgn.mins_);
                }
            } else if (last_action_.type_ == Action::Type::REMOTE_CHECK) {
                validator_.SendForwardResult(last_action_.forward_id_,
                        action_succeeded_);
            } else {
                // Simulation finished: we need to decide if we want to forward
                assert(action_succeeded_);
                validator_.adapter_->StopCollectingStats();
                const bool forwarded =
                        validator_.CheckForward(last_asgn_.get());
                if (!forwarded) {
                    last_action_.type_ = Action::Type::LOCAL_CHECK;
                    validator_.adapter_->SetAdapterMode(Adapter::EXACT);
                    return solver->RevAlloc(
                            new RestoreAssignment{last_asgn_.get()});
                }
            }
        }

        // check for searchlight termination
        if (validator_.CheckTerminate()) {
            // this will stop the search by failing the right branch
            LOG4CXX_INFO(logger, "Terminating the validator by force");
            solver->Fail();
        }

        // pick the next assignment
        if (asgns_.empty()) {
            // need to ask the validator
            Validator::CandidateVector *new_asgns =
                    validator_.GetNextAssignments();
            if (!new_asgns) {
                // No more assignments: stop the validator search
                LOG4CXX_INFO(logger, "Stopping the validator search");
                solver->Fail();
            }
            asgns_.swap(*new_asgns);

            LOG4CXX_TRACE(logger, "Got " << asgns_.size() << " new assignments "
                    "to check");
            delete new_asgns;
        }

        // Prepare assignment
        const int forw_id = asgns_.back().forw_id_;
        LiteToFullAssignment(*last_asgn_, asgns_.back().var_asgn_);
        asgns_.pop_back();

        last_action_.forward_id_ = forw_id;
        if (forw_id >= 0) {
            last_action_.type_ = Action::Type::REMOTE_CHECK;
        } else {
            last_action_.type_ = validator_.balancing_enabled_ ?
                    Action::Type::SIM : Action::Type::LOCAL_CHECK;
        }

        /*
         * Set adapter to either perform simulation or check immediately.
         */
        if (last_action_.type_ == Action::Type::SIM) {
            validator_.adapter_->SetAdapterMode(Adapter::DUMB);
            validator_.adapter_->StartCollectingStats();
        } else {
            validator_.adapter_->SetAdapterMode(Adapter::EXACT);
        }

        return solver->RevAlloc(new RestoreAssignment{last_asgn_.get()});
    }

    /**
     * Returns the name for debugging logging.
     * @return
     */
    virtual string DebugString() const {
        return "RestoreAssignment";
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
                my_fail_ = true;
                solver()->Fail();
            } else {
                if (!my_fail_) {
                    // Action wasn't successful -- fail
                    restore_db_.action_succeeded_ = false;
                }
                // Restart the search for garbage collecting
                my_fail_ = false;
                if (solver()->SearchDepth() > restart_period_) {
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

    // The validator producing the assignments
    Validator &validator_;

    // Assignments to validate
    Validator::CandidateVector asgns_;

    // True, if the last action was a success (no fail on solver)
    bool action_succeeded_;

    // The auxiliary monitor
    AuxRestoreMonitor aux_monitor_;

    // Assignment used in the last action
    AssignmentPtr last_asgn_;

    // Last action taken by the DB
    Action last_action_;
};

Validator::Validator(Searchlight &sl, SearchlightTask &sl_task,
        const StringVector &var_names) :
        sl_(sl),
        sl_task_(sl_task),
        solver_("validator solver"),
        adapter_(sl.CreateAdapter("validator")), // INTERVAL mode by default!
        search_vars_prototype_(&solver_) {

    // Start adapter in the dumb mode, since we don't need estimations
    adapter_->SetAdapterMode(Adapter::DUMB);

    // First, clone the solver
    if (!CloneModel(sl_, sl_.GetSolver(), solver_, adapter_)) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Cannot create the validator solver!";
    }

    // logging
    if (logger->isDebugEnabled()) {
        /*
         *  Unfortunately, it will dump it into std::cerr, which is
         *  redirected by SciDb to a file.
         */
        ModelVisitor *pmv = solver_.MakePrintModelVisitor();
        solver_.Accept(pmv);
    }

    // Create collector
    collector_ = solver_.RevAlloc(
            new SearchlightSolutionCollector(sl_task, &solver_));

    /*
     * Now, we want to find all decision variables. We do this via a
     * custom visitor.
     */
    VariableFinder var_finder;
    solver_.Accept(&var_finder);
    VariableFinder::StringVarMap var_map = var_finder.GetVarMap();
    for (StringVector::const_iterator cit = var_names.begin();
            cit != var_names.end(); cit++) {
        VariableFinder::StringVarMap::const_iterator var = var_map.find(*cit);
        if (var == var_map.end()) {
            std::ostringstream err_msg;
            err_msg << "Cannot find a variable in the duplicate solver: name="
                    << *cit;
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << err_msg.str();
        }
        /*
         *  const_cast is perfectly fine. Vars just come from the visitor,
         *  which does not change them. But these vars are supposed to be
         *  changed during the search anyway, albeit later.
         */
        search_vars_prototype_.Add(const_cast<IntVar *>(var->second));
        collector_->Add(const_cast<IntVar *>(var->second));
    }

    const SearchlightConfig &sl_config = sl.GetConfig();
    max_pending_validations_ =
            sl_config.get("searchlight.validator.max_validations", 1000000);
    restart_period_ = sl_config.get("searchlight.validator.restart_period",
            1024);

    balancing_enabled_ = sl_config.get("balance.validator_balance", 1);
    if (balancing_enabled_) {
        LOG4CXX_INFO(logger, "Forwarding is enabled for the validator...");
    }
}

void Validator::AddSolution(const Assignment &sol) {
    LiteVarAssignment lite_sol;
    FullAssignmentToLite(sol, lite_sol);

    LOG4CXX_TRACE(logger, "New solution to validate: " << sol.DebugString());

    std::unique_lock<std::mutex> validate_lock(to_validate_mtx_);
    to_validate_.push_back({std::move(lite_sol), -1});

    // We might have to wait -- flood control
    while (to_validate_.size() > max_pending_validations_) {
        /*
         * It is safe to use the same condition since the validator will be
         * certainly non-blocked.
         */
        LOG4CXX_TRACE(logger, "Waiting for the validator to catch up"
                ", queue size=" << to_validate_.size());
        validate_cond_.wait(validate_lock);
    }

    // Notify the validator in case it is blocked and unlock
    validate_cond_.notify_one();
}

void Validator::AddRemoteCandidates(LiteAssignmentVector &cands,
        InstanceID src, int forw_id) {
    std::lock_guard<std::mutex> validate_lock(to_validate_mtx_);
    while (!cands.empty()) {
        const int cand_id = remote_cand_id_++;
        to_validate_.push_back({std::move(cands.back()), cand_id});
        cands.pop_back();
        remote_candidates_.emplace(cand_id, std::make_pair(src, forw_id++));
    }
    validate_cond_.notify_one();
}

void Validator::HandleForwardResult(int id, bool result) {
    std::lock_guard<std::mutex> validate_lock{to_validate_mtx_};
    assert(forwarded_candidates_.find(id) != forwarded_candidates_.end());
    if (result) {
        sl_task_.ReportSolution(forwarded_candidates_[id].mins_);
    }
    forwarded_candidates_.erase(id);
    if (forwarded_candidates_.empty()) {
        validate_cond_.notify_one();
    }
}

void Validator::SendForwardResult(int forw_id, bool result) {
    std::unique_lock<std::mutex> validate_lock(to_validate_mtx_);
    const auto forw_info = remote_candidates_[forw_id];
    remote_candidates_.erase(forw_id);
    validate_lock.unlock();

    sl_task_.SendBalanceResult(forw_info.first, forw_info.second, result);
}

bool Validator::CheckForward(const Assignment *asgn) {
    const CoordinateSet &chunks = adapter_->GetCurrentStats().chunks_pos_;
    boost::shared_ptr<Query> query = sl_task_.GetQueryContext();
    std::vector<int> inst_counts(sl_task_.GetQueryInstanceCount(), 0);
    const InstanceID my_inst = sl_task_.GetInstanceID();

    // Get the distribution
    adapter_->GetSearchArrayDesc().GetChunksDistribution(
            query, chunks, inst_counts);

    /*
     *  For now we just pick the validator with the maximum count, which
     *  means it contains most of the chunks. In case the local instance
     *  contains the same number of chunks, we don't forward.
     */
    int max_inst = -1, max_count = -1;
    for (size_t i = 0; i < inst_counts.size(); ++i) {
        if (max_count < inst_counts[i]) {
            max_inst = i;
            max_count = inst_counts[i];
        }
    }
    if (max_inst != my_inst && inst_counts[my_inst] == max_count) {
        max_inst = my_inst;
    }

    if (max_inst != my_inst) {
        // forward
        LiteAssignmentVector asgns(1);
        FullAssignmentToLite(*asgn, asgns[0]);
        const int cand_id = forw_id_++;
        {
            std::lock_guard<std::mutex> lock{to_validate_mtx_};
            forwarded_candidates_.emplace(cand_id, asgns[0]);
        }
        sl_task_.ForwardCandidates(asgns, max_inst, cand_id);
        return true;
    }

    return false;
}

void Validator::Synchronize() const {
    std::unique_lock<std::mutex> validate_lock(to_validate_mtx_);
    while (!to_validate_.empty()) {
        LOG4CXX_TRACE(logger, "Synchronizing with the validator"
                ", queue size=" << to_validate_.size());
        validate_cond_.wait(validate_lock);
    }
}

void Validator::operator()() {
    LOG4CXX_INFO(logger, "Starting the validator search");
    DecisionBuilder *db =
            solver_.RevAlloc(new RestoreAssignmentBuilder(*this, &solver_,
                    restart_period_, &search_vars_prototype_));
    solver_status_ = solver_.Solve(db, collector_);
    LOG4CXX_INFO(logger, "Validator exited solve()");
}

Validator::CandidateVector *Validator::GetNextAssignments() {
    /*
     * For now, a simple policy: wait until we get an assignment to
     * validate and then check it.
     */
    std::unique_lock<std::mutex> validate_lock(to_validate_mtx_);
    while (true) {
        // Check interesting statuses first
        const Searchlight::Status sl_status = sl_.GetStatus();
        if (sl_status == Searchlight::Status::TERMINATED) {
            LOG4CXX_INFO(logger, "Terminating validator by force");
            // Notify the main solver about catching up, if needed
            to_validate_.clear();
            validate_cond_.notify_one();
            return nullptr;
        } else if (sl_status == Searchlight::Status::COMMITTED) {
            assert(FinishedLocally());
            LOG4CXX_INFO(logger, "Committing validator");
            return nullptr;
        } else if (sl_status == Searchlight::Status::FIN_SEARCH) {
            if (FinishedLocally()) {
                LOG4CXX_INFO(logger, "Validator finished locally");
                sl_.ReportFinValidator();
                // Still cannot exit: might have forwards coming up
            }
        }

        // Then check if we have any candidates or forwards
        if (!to_validate_.empty()) {
            break;
        }

        // Nothing new to do -- wait
        validate_cond_.wait(validate_lock);
    }

    // We should be here only if we have some candidates to check
    CandidateVector *res = new CandidateVector;
    res->swap(to_validate_);

    // Flow control: notify the main solver about catching up
    validate_cond_.notify_one();

    return res;
}

} /* namespace searchlight */

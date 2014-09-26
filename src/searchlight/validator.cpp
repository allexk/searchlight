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
     * @param validator servicing validator
     * @param adapter adapter used for access
     * @param s the solver
     * @param restart_period validator restart period
     * @param prototype prototype assignment to generate validations
     * @param asgns initial assignments to validate
     * @param slave true, if this is a slave builder for a helper
     */
    explicit RestoreAssignmentBuilder(Validator &validator,
            const AdapterPtr &adapter, Solver *s,
            int restart_period, const Assignment *prototype,
            Validator::CandidateVector &&asgns,
            bool slave)
        : validator_(validator),
          adapter_{adapter},
          asgns_{std::move(asgns)},
          action_succeeded_(false),
          aux_monitor_(*this, s, restart_period),
          last_asgn_{new Assignment{prototype}},
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
                adapter_->StopCollectingStats();
                const bool forwarded =
                        validator_.CheckForward(
                                adapter_->GetCurrentStats().chunks_pos_,
                                last_asgn_.get());
                if (!forwarded) {
                    last_action_.type_ = Action::Type::LOCAL_CHECK;
                    adapter_->SetAdapterMode(Adapter::EXACT);
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
            if (!slave_) {
                // need to ask the validator
                Validator::CandidateVector *new_asgns =
                        validator_.GetNextAssignments();
                if (!new_asgns) {
                    // No more assignments: stop the validator search
                    LOG4CXX_INFO(logger, "Stopping the validator search");
                    solver->Fail();
                }
                asgns_.swap(*new_asgns);

                LOG4CXX_TRACE(logger, "Got " << asgns_.size()
                        << " new assignments to check");
                delete new_asgns;
            } else {
                LOG4CXX_INFO(logger, "Helper validator finished workload...");
                solver->Fail();
            }
        }

        // Prepare assignment
        const int forw_id = asgns_.back().forw_id_;
        LiteToFullAssignment(*last_asgn_, asgns_.back().var_asgn_);
        asgns_.pop_back();

        last_action_.forward_id_ = forw_id;
        if (forw_id >= 0) {
            last_action_.type_ = Action::Type::REMOTE_CHECK;
        } else {
            last_action_.type_ =
                    validator_.forw_type_ != Validator::Forwarding::NONE ?
                    Action::Type::SIM : Action::Type::LOCAL_CHECK;
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

        return solver->RevAlloc(new RestoreAssignment{last_asgn_.get()});
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

    // True if this is a slave builder
    const bool slave_;
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
    }
    search_vars_prototype_.Store(); // Store initial values

    const SearchlightConfig &sl_config = sl_task_.GetConfig();
    max_pending_validations_ =
            sl_config.get("searchlight.validator.max_validations", 1000000);
    restart_period_ = sl_config.get("searchlight.validator.restart_period",
            1024);
    max_helpers_allowed_ =
            sl_config.get("searchlight.validator.max_helpers", 1);
    helper_workload_ =
            sl_config.get("searchlight.validator.helper_workload", 1000);

    // Determine forwarding type
    const std::string &forw_type =
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
}

void Validator::AddSolution(const Assignment &sol) {
    LiteVarAssignment lite_sol;
    FullAssignmentToLite(sol, lite_sol);

    LOG4CXX_TRACE(logger, "New solution to validate: " << sol.DebugString());

    std::unique_lock<std::mutex> validate_lock(to_validate_mtx_);
    PushCandidate({std::move(lite_sol), -1});

    // We might have to wait -- flood control
    while (to_validate_total_ > max_pending_validations_) {
        /*
         * It is safe to use the same condition since the validator will be
         * certainly non-blocked.
         */
        LOG4CXX_TRACE(logger, "Waiting for the validator to catch up"
                ", queue size=" << to_validate_total_);
        validate_cond_.wait(validate_lock);
    }

    // Notify the validator in case it is blocked and unlock
    validate_cond_.notify_one();
}

void Validator::PushCandidate(CandidateAssignment &&asgn) {
    if (to_validate_.empty() ||
            to_validate_.back().size() >= helper_workload_) {
        to_validate_.emplace_back();
        to_validate_.back().reserve(helper_workload_);
    }
    to_validate_.back().push_back(std::move(asgn));
    to_validate_total_++;
}

void Validator::AddRemoteCandidates(LiteAssignmentVector &cands,
        InstanceID src, int forw_id) {
    std::lock_guard<std::mutex> validate_lock(to_validate_mtx_);
    while (!cands.empty()) {
        const int cand_id = remote_cand_id_++;
        PushCandidate({std::move(cands.back()), cand_id});
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

bool Validator::CheckForward(const CoordinateSet &chunks,
        const Assignment *asgn) {
    const auto &active_validators = sl_task_.GetActiveValidators();
    std::vector<int> inst_counts(active_validators.size());

    // Get the distribution
    switch (forw_type_) {
        case Forwarding::STRIPES:
            adapter_->GetSearchArrayDesc().GetStripesChunkDistribution(
                    chunks, inst_counts);
            break;
        case Forwarding::DYNAMIC: {
            std::vector<int> inst_counts_all(sl_task_.GetQueryInstanceCount());
            adapter_->GetSearchArrayDesc().GetDynamicChunksDistribution(
                    sl_task_.GetQueryContext(), chunks, inst_counts_all);
            for (size_t i = 0; i < active_validators.size(); i++) {
                inst_counts[i] = inst_counts_all[active_validators[i]];
            }
            break;
        }
        case Forwarding::NONE:
            assert(false);
            break;
    }

    /*
     *  For now we just pick the validator with the maximum count, which
     *  means it contains most of the chunks.
     */
    int max_inst = -1, max_count = -1;
    for (size_t i = 0; i < inst_counts.size(); ++i) {
        if (max_count < inst_counts[i]) {
            max_inst = i;
            max_count = inst_counts[i];
        }
    }

    /*
     *  In case the local instance contains the same number of chunks,
     *  we don't forward.
     */
    if (my_logical_id_ != -1 && max_inst != my_logical_id_ &&
            inst_counts[my_logical_id_] == max_count) {

        max_inst = my_logical_id_;
    }

    if (max_inst != my_logical_id_) {
        // forward
        LiteAssignmentVector asgns(1);
        FullAssignmentToLite(*asgn, asgns[0]);
        const int cand_id = forw_id_++;
        {
            std::lock_guard<std::mutex> lock{to_validate_mtx_};
            forwarded_candidates_.emplace(cand_id, asgns[0]);
        }
        sl_task_.ForwardCandidates(asgns, active_validators[max_inst],
                cand_id);
        return true;
    }

    // This validator must be active here!
    assert(my_logical_id_ != -1);
    return false;
}

void Validator::Synchronize() const {
    std::unique_lock<std::mutex> validate_lock(to_validate_mtx_);
    while (!to_validate_.empty()) {
        LOG4CXX_TRACE(logger, "Synchronizing with the validator"
                ", queue size=" << to_validate_total_);
        validate_cond_.wait(validate_lock);
    }
}

void Validator::operator()() {
    LOG4CXX_INFO(logger, "Starting the validator search");
    DecisionBuilder *db =
            solver_.RevAlloc(new RestoreAssignmentBuilder(*this, adapter_,
                    &solver_, restart_period_, &search_vars_prototype_, {},
                    false));
    solver_status_ = solver_.Solve(db);
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

    // First, try to offload to helpers
    while (to_validate_.size() > 1) {
        ValidatorHelper *helper = DispatchValidatorHelper();
        if (helper) {
            CandidateVector &workload = to_validate_.front();
            LOG4CXX_INFO(logger, "Off-loading " << workload.size()
                    << " assignments to a helper...");
            to_validate_total_ -= workload.size();
            helper->RunWorkload(std::move(workload));
            to_validate_.pop_front();
        } else {
            break;
        }
    }

    // The rest we'll handle ourselves
    assert(to_validate_.size() > 0);
    CandidateVector *res = new CandidateVector;
    res->swap(to_validate_.front());
    to_validate_.pop_front();
    to_validate_total_ -= res->size();

    // Flow control: notify the main solver about catching up
    validate_cond_.notify_one();

    return res;
}

Validator::ValidatorHelper *Validator::DispatchValidatorHelper() {
    if (!free_validator_helpers_.empty()) {
        const int id = free_validator_helpers_.front();
        free_validator_helpers_.pop_front();
        ValidatorHelper *helper = validator_helpers_[id].get();
        return helper;
    } else if (validator_helpers_.size() < max_helpers_allowed_) {
        const int id = validator_helpers_.size();
        ValidatorHelper *helper =
                new ValidatorHelper{id, *this, search_vars_prototype_};
        validator_helpers_.emplace(id, ValidatorHelper::UniquePtr{helper});
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
                id_{id} {
    // Start adapter in the dumb mode, since we don't need estimations
    adapter_->SetAdapterMode(Adapter::DUMB);

    // First, clone the solver
    if (!CloneModel(parent_.sl_, parent_.solver_, solver_, adapter_)) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Cannot create a validator helper solver!";
    }

    /*
     * Now, we want to find all decision variables. We do this via a
     * custom visitor.
     */
    VariableFinder var_finder;
    solver_.Accept(&var_finder);
    VariableFinder::StringVarMap var_map = var_finder.GetVarMap();
    for (const auto &elem: init_assignment.IntVarContainer().elements()) {
        const std::string &var_name = elem.Var()->name();
        VariableFinder::StringVarMap::const_iterator it = var_map.find(var_name);
        if (it == var_map.end()) {
            std::ostringstream err_msg;
            err_msg << "Cannot find a variable in a duplicate solver: name="
                    << var_name;
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << err_msg.str();
        }
        /*
         *  const_cast is perfectly fine. Vars just come from the visitor,
         *  which does not change them. But these vars are supposed to be
         *  changed during the search anyway, albeit later.
         */
        IntVar *var = const_cast<IntVar *>(it->second);
        prototype_.Add(var);
        var->SetRange(elem.Var()->Min(), elem.Var()->Max());
    }
}

void Validator::ValidatorHelper::operator()() {
    LOG4CXX_INFO(logger, "Starting the helper validator search...");
    DecisionBuilder *db =
            solver_.RevAlloc(new RestoreAssignmentBuilder(parent_, adapter_,
                    &solver_, 0, &prototype_,
                    std::move(workload_), true));
    solver_.Solve(db);
    LOG4CXX_INFO(logger, "Validator helper exited solve()");

    std::lock_guard<std::mutex> lock{parent_.to_validate_mtx_};
    parent_.free_validator_helpers_.push_back(id_);
    parent_.validate_cond_.notify_one();
}

void Validator::ValidatorHelper::RunWorkload(CandidateVector &&workload) {
    assert(workload_.empty());
    Prepare();
    SetWorkload(std::move(workload));
    thr_ = std::thread{std::ref(*this)};
}

} /* namespace searchlight */

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
 * @file searchlight.cpp
 * The implementation of the main searchlight module.
 *
 * @author Alexander Kalinin
 */

#include "searchlight.h"
#include "validator.h"
#include "default_search.h"
#include "searchlight_task.h"

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.searchlight"));

void Searchlight::Cleanup() {
    // Wait for the solvers
    solvers_.clear();

    // Finally, destroy the validator
    if (validator_) {
        if (validator_thread_) {
            LOG4CXX_INFO(logger, "Waiting for the validator");
            validator_->WakeupIfIdle();
            validator_thread_->join();
            delete validator_thread_;
        }

        delete validator_;
    }
}

SearchlightSolver::~SearchlightSolver() {
    // Join the thread
    if (thr_.joinable()) {
        thr_.join();
    }

    // Time stats
    const auto secs = std::chrono::duration_cast<std::chrono::duration<double>>(
            total_solve_time_).count();
    LOG4CXX_INFO(logger, "Solver (0x" << std::hex << id_ << std::dec <<
            ") total time: " << secs << 's');

    /*
     * Search stats.
     *
     * Note that the solver's fail counter contains
     * not only true fails (due to violations), but also in-solution fails,
     * when we it just backtracks to find the next one. In our case, we always
     * continue after another candidate is found, so fails-sols should give
     * the number of true fails.
     */
    if (search_monitors_.validator_monitor_) {
        const int64 total_fails = solver_.failures();
        const int64_t total_candidates =
                search_monitors_.validator_monitor_->CandidatesNumber();
        LOG4CXX_INFO(logger, "Main search stats: fails=" << total_fails <<
                ", true fails=" << (total_fails - total_candidates) <<
                ", candidates=" << total_candidates);
    }

    if (status_ != Status::FINISHED) {
        LOG4CXX_INFO(logger, "Solver was terminated unexpectedly!");
    }
}

AttributeID Searchlight::RegisterAttribute(const std::string &name) {
    if (!array_desc_) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "No array registered with SL to register an attribute";
    }
    return array_desc_->RegisterAttribute(name,
            sl_task_.GetConfig().get("searchlight.load_aux_samples", 0));
}

const SearchlightConfig &Searchlight::GetConfig() const {
    return sl_task_.GetConfig();
}

void Searchlight::Prepare(const std::string &name, SLTaskFunc task_fun,
        uint64_t start_id, int solvers) {
    // Create solvers
    for (int i = 0; i < solvers; i++) {
        const uint64_t solver_id = start_id + i;
        solvers_.emplace_back(new SearchlightSolver(*this, solver_id, name));
        task_fun(this, i);
        // Note: task will call solver's prepare()
    }

    /*
     * Set up the validator. We start if the instance is active.
     * Depending on the active solvers/validators configuration, it
     * act as a full-fledged validator or just a forwarder.
     *
     * Even if the instance is in-active, we create it anyway, since it
     * simplifies message passing and maintenance (some messages are
     * broadcasted to be handled by everybody). Such a validator doesn't
     * consume any resources and pretty cheap to construct.
     *
     * TODO: We create the validator here, before the solver's Prepare()
     *  to catch the original model submitted by the user. There might be
     *  a better way to do this: export the model in a buffer and use the
     *  buffer for all subsequent model clones.
     */
    var_names_ = solvers_[0]->ComputeVarNames();
    LOG4CXX_INFO(logger, "Initiating the validator");
    validator_ = new Validator(*this, sl_task_, var_names_);
    if (sl_task_.InstanceActive(sl_task_.GetInstanceID())) {
        LOG4CXX_INFO(logger, "Starting the validator thread");
        validator_thread_ = new std::thread(std::ref(*validator_));

        // Prepare solvers
        for (auto &solver: solvers_) {
            solver->Prepare(*validator_);
        }
    } else {
        LOG4CXX_INFO(logger, "This instance is not active...");
        status_ = Status::COMMITTED;
    }
}

void SearchlightSolver::SetTask(const IntVarVector &primary_vars,
        const IntVarVector &secondary_vars,
        DecisionBuilder *db, const std::vector<SearchMonitor *> &monitors) {
    /*
     * Perform initialization
     */
    primary_vars_ = primary_vars;
    secondary_vars_ = secondary_vars;
    db_ = db;

    // All vars, for convenience
    all_vars_ = primary_vars_;
    all_vars_.insert(all_vars_.end(), secondary_vars_.begin(),
            secondary_vars_.end());

    // Store initial values to reset for future workloads
    vars_leaf_.Add(all_vars_);
    vars_leaf_.Store();

    // Save user monitors (the solver doesn't own them)
    search_monitors_.user_monitors_ = monitors;
}

void SearchlightSolver::Prepare(Validator &validator) {
    if (sl_.sl_task_.SolverActive(sl_.sl_task_.GetInstanceID())) {
        LOG4CXX_INFO(logger, "Solver is active at this instance...");
        // Establish monitors: validator (to transfer leaves) and terminator
        SearchLimit *terminator = solver_.MakeCustomLimit(
                NewPermanentCallback(this, &SearchlightSolver::CheckTerminate));
        search_monitors_.aux_monitors_.push_back(terminator);

        search_monitors_.validator_monitor_ = solver_.RevAlloc(
                new ValidatorMonitor{all_vars_, &solver_, validator});

        // Determine the worload; it will be assigned at Solve()
        DetermineLocalWorkload();

        // Switch adapters to interval mode
        for (auto &ad: adapters_) {
            ad->SetAdapterMode(Adapter::INTERVAL);
        }

        // logging model (only once)
        if (logger->isDebugEnabled()) {
            /*
             *  Unfortunately, it will dump it into std::cerr, which is
             *  redirected by SciDb to a file.
             */
            std::cerr << "Logging model for solver, id=0x" <<
                    std::hex << id_ << std::dec << "...\n";
            ModelVisitor *pmv = solver_.MakePrintModelVisitor();
            solver_.Accept(pmv);
        }

        // Decide if the balancing is enabled
        solver_balancing_enabled_ = sl_.GetConfig().get(
                "balance.solver_balance", 1);
        status_ = Status::PREPARED;
    } else {
        LOG4CXX_INFO(logger, "Solver is not configured at this instance...");
        status_ = Status::FINISHED;
    }
}

void SearchlightSolver::operator()() {
    try {
        // Solver always has some work to do if it's started
        bool work_expected = true;
        while (work_expected) {
            // Wait until something interesting happens
            {
                std::unique_lock<std::mutex> lock(mtx_);
                while (status_ == Status::VOID) {
                    wait_cond_.wait(lock);
                }
            }

            // Either solve or exit
            switch (status_) {
                case Status::PREPARED:
                    Solve();
                    if (status_ == Status::TERMINATED) {
                        LOG4CXX_DEBUG(logger, "Solver was terminated!");
                        work_expected = false;
                    }
                    break;
                case Status::FINISHED:
                case Status::TERMINATED:
                    work_expected = false;
                    break;
                default:
                    assert(false);
                    throw SYSTEM_EXCEPTION(
                            SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                            << "Unexpected status in the solver loop!";
                    break;
            }
        }
    } catch (const scidb::Exception &ex) {
        /*
         * std::thread would call std::terminate(), so
         *  we have to take care of proper error reporting
         */
        sl_.sl_task_.HandleSearchlightError(ex.copy());
    } catch (const std::exception &e) {
        // Catch other C++ and library exceptions and translate them
        sl_.sl_task_.HandleSearchlightError((SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                SCIDB_LE_ILLEGAL_OPERATION) << e.what()).copy());
    } catch (...) {
        sl_.sl_task_.HandleSearchlightError((SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                SCIDB_LE_ILLEGAL_OPERATION) <<
                "Unknown exception in SL!").copy());
    }
}

StringVector SearchlightSolver::ComputeVarNames() const {
    // All vars, for convenience
    StringVector var_names(all_vars_.size());
    int i = 0;
    for (IntVarVector::const_iterator cit = all_vars_.begin();
            cit != all_vars_.end(); cit++) {
        const IntVar * const var = *cit;
        std::string var_name;
        if (var->HasName()) {
            var_name = var->name();
        } else {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                    SCIDB_LE_ILLEGAL_OPERATION)
                    << "Every decision variable must have a name!";
        }
        var_names[i++] = var_name;
    }

    return var_names;
}

void SearchlightSolver::PrepareHelper(LiteAssignmentVector &load) {
    std::lock_guard<std::mutex> lock{mtx_};

    assert(status_ == Status::VOID);
    assert(helper_load_.empty());

    LOG4CXX_INFO(logger, "Preparing solver for remote work...");
    helper_load_.swap(load);

    // Decide if the balancing is enabled
    solver_balancing_enabled_ =
            sl_.sl_task_.GetConfig().get("balance.solver_balance", 1);
    status_ = Status::PREPARED;

    wait_cond_.notify_one();
}

void SearchlightSolver::DetermineLocalWorkload() {
    /*
     * Determine our partition. First, determine the variable with the
     * largest domain and the assume equivalent stripes along its domain.
     */
    IntVar *split_var = nullptr;
    uint64_t max_size = 0;
    size_t split_var_num = 0;
    for (size_t i = 0; i < primary_vars_.size(); i++) {
        const auto var = primary_vars_[i];
        const uint64_t var_size = var->Max() - var->Min() + 1;
        if (!split_var || var_size > max_size) {
            split_var = var;
            max_size = var_size;
            split_var_num = i;
        }
    }
    assert(split_var);

    /*
     * Then, determine intervals.
     */
    const size_t active_solvers_num = sl_.sl_task_.GetTotalSolversNum();
    const uint32_t solver_id = sl_.sl_task_.GetGlobalOrdinalSolverId(id_);

    int slices_number =
            sl_.sl_task_.GetConfig().get("balance.slices_number", 0);
    if (slices_number < active_solvers_num || active_solvers_num == 1) {
        slices_number = active_solvers_num;
    } else if (slices_number > split_var->Size()) {
        slices_number = split_var->Size();
    }
    const uint64_t slice_len = ceil(double(max_size) / slices_number);
    LOG4CXX_INFO(logger, "The number of slices is " << slices_number
            << " and the slice length is " << slice_len);
    const bool load_slices =
            sl_.sl_task_.GetConfig().get("balance.load_slices", 1);
    if (load_slices) {
        // See comment in Solve() about NewSearch()
        solver_.NewSearch(db_, search_monitors_.GetSearchMonitors());
        vars_leaf_.FreezeQueue();
        // Remove first hole
        const int64_t initial_var_min = split_var->Min();
        split_var->RemoveInterval(split_var->Min(),
                split_var->Min() + slice_len * solver_id - 1);

        // Remaining holes (note, Min() might have changed!)
        int64_t hole_start = initial_var_min + slice_len * (solver_id + 1);
        while (hole_start <= split_var->Max()) {
            const int64_t hole_end =
                    hole_start + (active_solvers_num - 1) * slice_len - 1;
            split_var->RemoveInterval(hole_start, hole_end);
            hole_start = hole_end + 1 + slice_len;
        }
        vars_leaf_.UnfreezeQueue();
        LOG4CXX_INFO(logger, "Loaded slices into the split var: "
                << split_var->DebugString());
        local_remote_override_ = true;
    } else {
        for (int64_t start = split_var->Min() + slice_len * solver_id;
                start <= split_var->Max();
                start += active_solvers_num * slice_len) {
            // end might fall out of range, which is okay
            const int64_t end = start + slice_len - 1;

            // we define the load in terms of initial vars + split_var intervals
            LiteVarAssignment asgn;
            FullAssignmentToLite(vars_leaf_, asgn);
            asgn.mins_[split_var_num] = start;
            asgn.maxs_[split_var_num] = end;
            local_load_.push_back(std::move(asgn));
        }
    }
}

void SearchlightSolver::Solve() {
    // Enter the search
    status_ = Status::SEARCHING;

    // Initially we always have work: local or remote
    bool solver_has_work = !local_load_.empty() || !helper_load_.empty() ||
            local_remote_override_;
    while (solver_has_work) {
        if (!local_remote_override_) {
            // Local load has priority
            if (!local_load_.empty()) {
                LiteToFullAssignment(vars_leaf_, local_load_.back());
                local_load_.pop_back();

            } else if (!helper_load_.empty()) {
                // Get a remote assignment from the remote load
                LiteToFullAssignment(vars_leaf_, helper_load_.back());
                helper_load_.pop_back();
            } else {
                assert(false);
            }

            /*
             * Here we want to create a new search. It is crucial to be able to
             * roll-back all changes when this solver gets a new local or remote
             * assignment. The only way to do this is to
             * initialize the search now, so that the Solver puts a
             * INITIAL_SEARCH marker in the trace. Then, after the Solve() is
             * finished it will roll-back everything up to this point.
             */
            solver_.NewSearch(db_, search_monitors_.GetSearchMonitors());
            LOG4CXX_INFO(logger, "Setting up solver assignment: "
                    << vars_leaf_.DebugString());
            vars_leaf_.Restore();
        } else {
            local_remote_override_ = false;
        }

        // Starting the timer
        const auto solve_start_time = std::chrono::steady_clock::now();

        /*
         * Since the collector will visit all the leaves, the loop is empty.
         * Also, some heuristics, like SLRandom, find candidates via
         * nested search, which is not going to show up here anyway.
         */
        LOG4CXX_INFO(logger, "Starting the main search");
        while (solver_.NextSolution());
        solver_.EndSearch();

        // stopping the timer
        const auto solve_end_time = std::chrono::steady_clock::now();
        total_solve_time_ +=
                std::chrono::duration_cast<decltype(total_solve_time_)>(
                solve_end_time - solve_start_time);

        // Continue only if we have another assignment to set up
        solver_has_work = !local_load_.empty() || !helper_load_.empty() ||
                local_remote_override_;
    }

    LOG4CXX_INFO(logger, "Finished solver's workload...");
    if (status_ != Status::TERMINATED) {
        status_ = Status::VOID;

        // Report idleness
        sl_.sl_task_.ReportIdleSolver(id_, true);

        // Then get rid of remaining helpers
        RejectHelpers(false);
    }

    /*
     *  If we're terminated, we don't have to clean up:
     *  the system will go down anyway without blocking.
     */
}

bool ValidatorMonitor::AtSolution() {
    Assignment * const asgn = prototype_.get();

    // Store the solution (assuming complete assignment)
    asgn->Store();
    LOG4CXX_TRACE(logger, "Encountered a leaf: " << asgn->DebugString() <<
            ", depth=" << solver()->SearchDepth());

    // TODO: check if need this; a leaf should always be complete
    bool complete = true;
    for (const auto &var: asgn->IntVarContainer().elements()) {
        if (!var.Bound()) {
            complete = false;
            break;
        }
    }

    assert(complete);
    if (complete) {
        candidates_++;
        validator_.AddSolution(*asgn);
    }

    return true;
}

void Searchlight::ReportFinValidator() {
    status_ = Status::FIN_VALID;
    sl_task_.ReportFinValidator();
}

void Searchlight::HandleEndOfSearch() {
    status_ = Status::FIN_SEARCH;

    // Solvers
    for (auto &s: solvers_) {
        s->HandleEndOfSearch();
    }

    // Validator
    assert(validator_);
    validator_->WakeupIfIdle();
}

void Searchlight::HandleCommit() {
    status_ = Status::COMMITTED;
    assert(validator_);
    validator_->WakeupIfIdle();
}

void SearchlightSolver::HandleHelper(uint64_t id) {
    {
        std::lock_guard<std::mutex> lock(mtx_);
        helpers_.push_back(id);
    }

    // Check for reject
    if (status_ != Status::SEARCHING && status_ != Status::PREPARED) {
        RejectHelpers(false);
    } else if (!solver_balancing_enabled_) {
        RejectHelpers(true);
    }
}

void SearchlightSolver::RejectHelpers(bool hard) {
    std::unique_lock<std::mutex> lock(mtx_);
    std::vector<uint64_t> inst(helpers_.begin(), helpers_.end());
    helpers_.clear();
    lock.unlock();

    if (!inst.empty()) {
        sl_.sl_task_.RejectHelp(inst, id_, hard);
    }
}

void SearchlightSolver::DispatchWork(
        const LiteAssignmentVector &work) {
    assert(HelpAvailable());
    InstanceID helper;
    {
        std::lock_guard<std::mutex> lock(mtx_);
        helper = helpers_.front();
        helpers_.pop_front();
    }
    sl_.sl_task_.DispatchWork(work, helper);
}

std::string Searchlight::SolutionToString(
        const std::vector<int64_t> &vals) const {
    // stringify the solution
    std::ostringstream sol_string;
    for (size_t i = 0; i < var_names_.size(); i++) {
        sol_string << var_names_[i] << "=" << vals[i];
        if (i != var_names_.size() - 1) {
            sol_string << ", ";
        }
    }
    return sol_string.str();
}

DecisionBuilder *SearchlightSolver::CreateDefaultHeuristic(
        const IntVarVector &primary_vars,
        const IntVarVector &secondary_vars) {
    return solver_.RevAlloc(
            new SLSearch(sl_, *this, primary_vars, secondary_vars));
}

SearchMonitor *SearchlightSolver::CreateBalancingMonitor(
        const IntVarVector &vars,
        double low, double high) {
    LOG4CXX_INFO(logger, "Creating general balancer with interval ["
            << low << ", " << high << "]");
    return solver_.RevAlloc(
            new BalancingMonitor(*this, sl_, vars, low, high));
}

AdapterPtr SearchlightSolver::CreateAdapter(const std::string &name) {
    auto res = sl_.CreateAdapter(name + "_" + std::to_string(id_));
    adapters_.push_back(res);
    return res;
}

SearchMonitor *MakeCumulativeTimeLimit(Solver &s, int64 time_ms) {
    return s.MakeLimit(time_ms, kint64max, kint64max, kint64max, true, true);
}

} /* namespace searchlight */

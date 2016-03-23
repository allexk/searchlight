/* Copyright 2014, Brown University, Providence, RI.
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
 * @file default_search.cpp
 *
 * The implementation of the default search.
 *
 * @author Alexander Kalinin
 */

#include "default_search.h"
#include "validator.h"

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.sl_heuristic"));

/**
 * Prints the impact to the output stream.
 *
 * @param s the stream
 * @param impact the impact
 * @return the same stream
 */
std::ostream &operator<<(std::ostream &s, const SLSearch::Impact &impact) {
    s << "(samples: " << impact.tried_assigns_;
    s << ", succ.: " << impact.succ_assigns_;
    s << ", chunks: " << impact.chunks_to_read_;
    s << ", penalty: " << impact.penalty_;
    s << ", active: " << impact.active;
    s << ')';

    return s;
}

namespace {

/**
 * This monitor tracks the progress of the interval explorer and collects
 * relevant statistics.
 *
 * This monitor also serves as a limiter that finishes the exploring search
 * after the specified number of fails.
 */
class IntervalImpactMonitor : public SearchMonitor {
public:
    /**
     * Creates a new explorer monitor.
     *
     * @param s the solver for the search
     */
    IntervalImpactMonitor(Solver * const s) :
        SearchMonitor(s),
        fails_(0),
        leaves_(0),
        fails_limit_(0) {}

    /**
     * Set the limit of fails after which the monitor will fail the search.
     *
     * @param limit fails limit
     */
    void SetFailsLimit(int limit) {
        fails_limit_ = limit;
    }

    /**
     * Returns the number of true fails. A true fail is any fail except
     * a backtracking from a leaf.
     *
     * @return the number of true fails
     */
    int GetTrueFails() const {
        assert(fails_ >= leaves_);
        return fails_ - leaves_;
    }

    /**
     * Returns the total number of fails.
     *
     * @return the total number of fails
     */
    int GetTotalFails() const {
        return fails_;
    }

    /**
     * Returns the number of leaves (solutions) reached.
     *
     * @return the number of leaves reached
     */
    int GetLeavesReached() const {
        return leaves_;
    }

    /**
     * Callback called when a fail begins, before the jump.
     */
    virtual void BeginFail() override {
        fails_++;
        if (fails_ >= fails_limit_) {
            FinishCurrentSearch();
            // Not calling Fail() since we're failing anyway
        }
    }

    /**
     * Callback called when a leaf is reached.
     *
     * @return true, if we want to backtrack and continue; false, otherwise
     */
    virtual bool AtSolution() override {
        leaves_++;
        return true;
    }

private:
    // Total fails
    int fails_;

    // Leaves reached
    int leaves_;

    // The maximum number of fails
    int fails_limit_;
};

/**
 * This monitors checks the number of fails and leaves reached periodically,
 * and finishes the search if the ratio is below the threshold.
 */
class FinishOnFailsMonitor : public SearchMonitor {
public:
    /**
     * Creates a new check-fail-and-finish monitor.
     *
     * @param s the solver
     * @param period the length of checking period in fails
     * @param threshold the threshold of finish
     */
    FinishOnFailsMonitor(Solver *s, int64_t period, double threshold) :
        SearchMonitor(s),
        fails_period_(period), finish_threshold_(threshold) {}

    /**
     * Called when the solver begins the fail.
     */
    virtual void BeginFail() override {
        fails_++;
        if (fails_ == fails_period_) {
            const double ratio = double(leaves_) / fails_;
            if (ratio <= finish_threshold_) {
                FinishCurrentSearch();
                LOG4CXX_TRACE(logger, "Finishing the search because "
                        "of a large fail ratio: " << ratio);
            }

            fails_ = leaves_ = 0;
        }
    }

    /**
     * Called at a leaf.
     *
     * @return true, since we want to continue
     */
    virtual bool AtSolution() override {
        leaves_++;
        return true;
    }

private:
    // The length of checking period in fails
    const int64_t fails_period_;

    // Fails since the last check
    int64_t fails_ = 0;

    // The number of leaves since the last check
    int64_t leaves_ = 0;

    // Threshold of leaves/fails at which the search will be finished
    const double finish_threshold_;
};

}

/**
 * This class represents a decision to set an interval for a variable.
 */
class SLSearch::SetIntervalDecision : public Decision {
public:
    /**
     * Creates a new decision that sets the interval.
     *
     * @param search the SL search builder
     * @param var_int_ind the interval's index in the builder
     * @param var the variable to set
     * @param min the left interval bound
     * @param max the right interval bound
     */
    SetIntervalDecision(SLSearch &search, int var_int_ind, IntVar *var,
            int64 min, int64 max) :
                sl_search_(search), var_int_ind_(var_int_ind),
                var_(var), min_(min), max_(max) {}

    /**
     * The left branch. In this case it sets the interval.
     *
     * @param s the solver
     */
    virtual void Apply(Solver* const s) override {
        LOG4CXX_TRACE(logger, "Setting interval for  " <<
                var_->DebugString() << " to: [" << min_ << ", " << max_ << "]");
        var_->SetRange(min_, max_);
        LOG4CXX_TRACE(logger, "Interval set");
    }

    /**
     * The right branch. It removes the interval.
     *
     * @param s the solver
     */
    virtual void Refute(Solver* const s) override {
        LOG4CXX_TRACE(logger, "Temporarily deactivating interval for  " <<
                var_->DebugString() << ", [" << min_ << ", " << max_ << "]");
        /*
         * Temporarily freeze the solver's queue, since some variables
         * perform element-wise removes even for intervals.
         */
        var_->FreezeQueue();
        var_->RemoveInterval(min_, max_);
        var_->UnfreezeQueue();
        // Temporarily deactivate the interval
        sl_search_.DeactivateVarInterval(var_, var_int_ind_, true);
        LOG4CXX_TRACE(logger, "Interval removed and deacrivated");
    }

    /**
     * Returns a string for debug logging.
     *
     * @return a debug string
     */
    virtual std::string DebugString() const  override {
        return "SetIntervalDecision";
    }

private:
    // The SL heuristic
    SLSearch &sl_search_;

    // The interval index of the variable
    const int var_int_ind_;

    // The variable to set
    IntVar * const var_;

    // The left interval bound
    const int64 min_;

    // The right interval bound
    const int64 max_;
};

/**
 * Search heuristic to randomly sample the specified interval for the
 * specified number of tries (fails).
 */
class SLSearch::IntervalImpactBuilder : public DecisionBuilder {
public:
    /**
     * Creates a new interval sample heuristic.
     *
     * @param sl_search top-level SL search
     * @param impact_monitor impact monitor that tracks fails
     * @param var the var to bind to an interval
     * @param min the left interval bound
     * @param max the right interval bound
     * @param monitors monitors to establish, vector will be modified inside,
     *     but left unchanged at exit
     */
    IntervalImpactBuilder(const SLSearch &sl_search,
            IntervalImpactMonitor &impact_monitor,
            IntVar * const var,
            const int64_t min, const int64_t max,
            std::vector<SearchMonitor *> &monitors) :
            sl_search_(sl_search),
            impact_monitor_(impact_monitor),
            var_(var), min_(min), max_(max),
            monitors_(monitors),
            problem_infeasible_(true) {}

    /**
     * Produces a new decision. In this case it sets the interval and initiates
     * a nested search to sample it.
     *
     * @param s the solver
     * @return the next decision; nullptr is at a leaf
     */
    virtual Decision* Next(Solver* const s) override {
        LOG4CXX_TRACE(logger, "Sampling interval for " << var_->DebugString()
                << ", [" << min_ << ", " << max_ << "]");

        var_->SetRange(min_, max_);
        problem_infeasible_ = false;

        /*
         *  Determine number of probes.
         *
         *  Default: percentage of the size of the corresponding
         *  hyper-rectangle with fixed interval for this variable.
         *
         *  Next: we limit it by the maximum number (from config).
         */
        uint64_t ss_size = 1;
        for (const IntVar *var: sl_search_.all_vars_) {
            ss_size *= var->Size();
        }
        int probes_limit = ss_size *
                sl_search_.search_config_.probes_percentage_;
        probes_limit = std::min(probes_limit,
                sl_search_.search_config_.max_probes_number_);
        impact_monitor_.SetFailsLimit(probes_limit);

        // We are going to conduct a random search...
        DecisionBuilder * const random_db = s->MakePhase(sl_search_.all_vars_,
                Solver::CHOOSE_RANDOM, Solver::ASSIGN_RANDOM_VALUE);

        // Nested search
        int luby_scale = sl_search_.search_config_.luby_scale_;
        if (luby_scale != 0) {
            LOG4CXX_TRACE(logger, "Setting Luby restarts for impact estimator");
            monitors_.push_back(s->MakeLubyRestart(luby_scale));
        }
        s->Solve(random_db, monitors_);
        if (luby_scale != 0) {
            monitors_.pop_back();
        }

        // After the nested search, we have nothing to do. Exit
        return nullptr;
    }

    /**
     * Determines if the problem was infeasible.
     *
     * @return true, if the problem was infeasible; false, otherwise
     */
    bool ProblemInfeasible() const {
        return problem_infeasible_;
    }

    /**
     * Returns a debug string.
     *
     * @return string for debugging
     */
    virtual std::string DebugString() const override {
        return "Searchlight Interval Explorer Search";
    }

private:
    // Top-level SL search
    const SLSearch &sl_search_;

    // Impact monitor that tracks fails
    IntervalImpactMonitor &impact_monitor_;

    // The variable to set
    IntVar * const var_;

    // The interval to set
    const int64_t min_, max_;

    // Monitors to establish for the explorer search
    std::vector<SearchMonitor *> &monitors_;

    // True, if the problem was infeasible (i.e., the interval is infeasible)
    bool problem_infeasible_;
};

SLSearch::SLSearch(Searchlight &sl,
        SearchlightSolver &sl_solver,
        const IntVarVector &primary_vars,
        const IntVarVector &secondary_vars) :
    primary_vars_(primary_vars),
    secondary_vars_(secondary_vars),
    intervals_explored_(false),
    sl_(sl),
    sl_solver_(sl_solver),
    solver_(sl_solver.GetSearchSolver()),
    dummy_monitor_(&solver_),
    search_config_(sl.GetConfig()) {

    // Init impacts
    for (auto int_var: primary_vars_) {
        if (!int_var->Bound()) {
            const int var_ind = var_to_index_.size();
            var_to_index_[int_var] = var_ind;

            const int64 var_min = int_var->Min();
            const int64 var_max = int_var->Max();
            const int64 var_len = var_max - var_min + 1;
            const size_t interval_len = var_len <= search_config_.splits_ ?
                    1 : var_len / search_config_.splits_;
            var_impacts_.push_back(IntervalVarImpacts(var_min, var_max,
                    interval_len, int_var));
        }
    }

    all_vars_.reserve(primary_vars.size() + secondary_vars.size());
    all_vars_.insert(all_vars_.end(), primary_vars.begin(),
            primary_vars.end());
    all_vars_.insert(all_vars_.end(), secondary_vars.begin(),
            secondary_vars.end());

    if (logger->isDebugEnabled()) {
        std::ostringstream deb;
        deb << "SLSearch config follows:\n";
        OutputConfig(deb);
        logger->debug(deb.str(), LOG4CXX_LOCATION);
    }
}

Decision* SLSearch::Next(Solver* const s) {
    if (!intervals_explored_ && !VarsIntervalBound()) {
        const auto init_start_time = std::chrono::steady_clock::now();
        InitIntervals(s);
        s->SaveAndSetValue(&intervals_explored_, true);
        //intervals_explored_ = true;
        const auto init_end_time = std::chrono::steady_clock::now();
        const int64_t init_seconds =
                std::chrono::duration_cast<std::chrono::seconds>
                (init_end_time - init_start_time).count();
        LOG4CXX_INFO(logger, "The init went for " << init_seconds << "s.");
    }

    // Here, we have all impacts computed
    if (!VarsIntervalBound()) {
        const VarImpactInfo best_var_imp = FindBestVar();
        if (best_var_imp.var_index >= 0) {
            const int var_ind = best_var_imp.var_index;
            const int int_ind = best_var_imp.int_index;

            IntVar * const var = var_impacts_[var_ind].var_;
            const Impact &impact = var_impacts_[var_ind].impacts_[int_ind].second;
            const int64 int_min = var_impacts_[var_ind].impacts_[int_ind].first;
            const int64 int_max = int_min +
                    var_impacts_[var_ind].interval_length_ - 1;

            LOG4CXX_TRACE(logger, "Determined the next best interval for " <<
                    var << ", [" << int_min << ", " <<
                    int_max << "]" << ", Impact: " << impact);

            // we should penalize the interval
            var_impacts_[var_ind].impacts_[int_ind].second.penalty_ *=
                    search_config_.interval_penalty_;

            // and try to set it
            return s->RevAlloc(new SetIntervalDecision(*this, int_ind,
                    var, int_min, int_max));
        } else {
            // cannot find a suitable interval: exhausted
            LOG4CXX_TRACE(logger,
                    "Could not find the next best interval: exhausted");
            return s->MakeFailDecision();
        }
    }

    /*
     *  Check if we can off-load some work. We can if:
     *    1) We are not exploring the last region
     *    2) We have some help available
     *
     *  For now, we just send out a single region.
     */
    if (s->SearchLeftDepth() > 0 && sl_solver_.HelpAvailable()) {
        // Take a snapshot of vars
        AssignmentPtr asgn{new Assignment(s)};
        asgn->Add(all_vars_);
        asgn->Store();

        // Log
        LOG4CXX_DEBUG(logger, "Off-loading a region to a helper: "
                << asgn->DebugString());

        CandidateVector work(1);
        FullAssignmentToLite(*asgn, work.back().var_asgn_);
        sl_solver_.DispatchWork(work);

        return s->MakeFailDecision();
    }

    LOG4CXX_DEBUG(logger, "Initiating an in-interval search");

    // All variables are bound to intervals -- random search with restarts
    // We are going to conduct a random search...
    DecisionBuilder * const random_db = s->MakePhase(all_vars_,
            Solver::CHOOSE_RANDOM, Solver::ASSIGN_RANDOM_VALUE);

    // set the time limit if needed
    std::vector<SearchMonitor *> monitors(sl_solver_.GetAuxMonitors());
    const auto &user_monitors = sl_solver_.GetUserMonitors();
    monitors.insert(monitors.end(), user_monitors.begin(),
            user_monitors.end());
    monitors.push_back(sl_solver_.GetValidatorMonitor());

    if (search_config_.time_strategy_ == SLConfig::CONST) {
        monitors.push_back(s->MakeTimeLimit(
                search_config_.time_interval_ * 1000));
        /*
         * Make Luby restarts for better coverage.
         *
         * Note, we don't establish Luby restarts if we don't have a time limit,
         * since the search will run infinitely long.
         */
        if (search_config_.luby_scale_ != 0) {
            monitors.push_back(s->MakeLubyRestart(search_config_.luby_scale_));
        }
    }

    /*
     * Establish a finish-on-frequent fails monitor. The rationale behind
     * this is that if we fail often before a leaf is reached, we are
     * probably wasting our time, since the search is random. If we don't
     * have a time limit, we might finish this search faster.
     *
     * Note, we don't establish Luby restarts if we don't have a time limit,
     * since the search will run infinitely long.
     */
    monitors.push_back(s->RevAlloc(new FinishOnFailsMonitor{s,
        search_config_.fails_restart_probes_,
        search_config_.fails_restart_thr_}));

    // starting the timer
    const auto solve_start_time = std::chrono::steady_clock::now();
    s->Solve(random_db, monitors);
    const auto solve_end_time = std::chrono::steady_clock::now();

    // time elapsed (might be less than a half of search_time_limit_, if
    // FinishOnFails finished the search.
    int64_t solve_seconds =
            std::chrono::duration_cast<std::chrono::seconds>
            (solve_end_time - solve_start_time).count();
    LOG4CXX_DEBUG(logger, "The search went for " << solve_seconds << "s.");

    // wait for the validator
    if (search_config_.validator_synchronize_) {
        LOG4CXX_DEBUG(logger, "Waiting for the validator...");
        sl_.GetValidator().Synchronize();
    }

    /*
     *  Will try a different interval combination. Since we penalize
     *  intervals we choose we should end up with a different combination,
     *  unless the current intervals are _reall_ promising (having a
     *  high impact).
     */
    if (search_config_.do_restarts_) {
        LOG4CXX_DEBUG(logger, "In-interval search, restarting.");
        dummy_monitor_.RestartCurrentSearch();
    }

    return s->MakeFailDecision();
}

void SLSearch::InitIntervals(Solver * const s) {
    LOG4CXX_INFO(logger, "Exploring interval impacts");
    std::vector<SearchMonitor *> nested_monitors(sl_solver_.GetAuxMonitors());
    if (search_config_.submit_probes_) {
        nested_monitors.push_back(sl_solver_.GetValidatorMonitor());
    }

    for (auto &var_impact: var_impacts_) {
        if (var_impact.impacts_.size() > 1) { // multiple intervals
            for (auto &interval_impact: var_impact.impacts_) {
                Impact &imp = interval_impact.second;
                const int64 min = interval_impact.first;
                const int64 max = min + var_impact.interval_length_ - 1;

                // Does the interval belong to the domain?
                const int64 var_min = var_impact.var_->Min();
                const int64 var_max = var_impact.var_->Max();
                if (var_min > max || var_max < min) {
                    LOG4CXX_TRACE(logger, "The interval is outside of domain");
                    s->SaveAndSetValue(&imp.active, false);
                    continue;
                }

                LOG4CXX_TRACE(logger, "Exploring interval [" <<
                        min << ", " << max << "] for " <<
                        var_impact.var_->DebugString());

                IntervalImpactMonitor explorer_monitor{s};
                nested_monitors.push_back(&explorer_monitor);
                IntervalImpactBuilder explorer{*this, explorer_monitor,
                    var_impact.var_, min, max, nested_monitors};
                s->Solve(&explorer);
                nested_monitors.pop_back();

                s->SaveAndSetValue(&imp.tried_assigns_,
                        explorer_monitor.GetTotalFails());
                s->SaveAndSetValue(&imp.succ_assigns_,
                        explorer_monitor.GetLeavesReached());
                if (explorer.ProblemInfeasible()) {
                    // We should deactivate the interval
                    LOG4CXX_TRACE(logger, "The interval is infeasible");
                    s->SaveAndSetValue(&imp.active, false);
                }
                LOG4CXX_TRACE(logger, "Determined impact for " <<
                        "[" << min << ", " << max << "]: " << imp);
            }
        }
    }
    LOG4CXX_INFO(logger, "Finished exploring interval impacts");
    if (logger->isDebugEnabled()) {
        std::ostringstream deb;
        deb << "SLSearch impacts for primary variables follow:\n";
        OutputImpacts(deb);
        logger->debug(deb.str(), LOG4CXX_LOCATION);
    }

    // wait for the validator
    if (search_config_.validator_synchronize_) {
        LOG4CXX_DEBUG(logger, "Waiting for the validator...");
        sl_.GetValidator().Synchronize();
    }
}

SLSearch::VarImpactInfo SLSearch::FindBestVar() const {
    VarImpactInfo res;

    // Choose from all variables
    for (size_t var_ind = 0; var_ind < var_impacts_.size(); var_ind++) {
        const IntervalVarImpacts &var_intervals = var_impacts_[var_ind];
        const IntVar * const var = var_intervals.var_;

        const int64_t var_min = var->Min();
        const int64_t var_max = var->Max();
        const int64_t domain_len = var_max - var_min + 1;
        const int64_t var_int_len = var_intervals.interval_length_;

        if (domain_len > var_int_len) {
            // The domain is large -- the var should be bound to an interval
            double best_impact = -2;
            int best_int_ind = -1;

            // Choose from all intervals
            for (size_t int_ind = 0; int_ind < var_intervals.impacts_.size();
                    int_ind++) {
                const IntervalImpact &interval =
                        var_intervals.impacts_[int_ind];
                // Do the domain range and the interval intersect?
                const int64 int_min = interval.first;
                const int64 int_max = interval.first + var_int_len - 1;
                if (int_max >= var_min && int_min <= var_max &&
                        interval.second.active) {
                    // Intervals intersect and the impact is active
                    // TODO: try other types of impact
                    const double imp = interval.second.PromiseImpact();
                    if (imp > best_impact) {
                        best_impact = imp;
                        best_int_ind = int_ind;
                    }
                }
            }
            if (best_impact > res.impact_) {
                res.var_index = var_ind;
                res.int_index = best_int_ind;
                res.impact_ = best_impact;
            }
        }
    }

    return res;
}

bool SLSearch::VarsIntervalBound() const {
    bool bound = true;
    for (const auto &var_impact: var_impacts_) {
        const IntVar * const var = var_impact.var_;
        const int64_t domain_len = var->Max() - var->Min() + 1;
        const int64_t var_int_len = var_impact.interval_length_;

        if (domain_len > var_int_len) {
            bound = false;
            break;
        }
    }
    return bound;
}

std::ostream &SLSearch::OutputImpacts(std::ostream &os) const {
    for (const auto &var_imps: var_impacts_) {
        const size_t len = var_imps.interval_length_;
        os << "Impacts for " << var_imps.var_->DebugString() << '\n';
        for (const auto &imp: var_imps.impacts_) {
            const int64_t min = imp.first;
            const int64_t max = min + len - 1;
            os << "\t[" << min << ", " << max << "]: " << imp.second << '\n';
        }
    }

    return os;
}

std::ostream &SLSearch::OutputConfig(std::ostream &os) const {
    os << "Time strategy: ";
    switch (search_config_.time_strategy_) {
        case SLConfig::CONST:
            os << "constant: " << search_config_.time_interval_ << " seconds";
            break;
        case SLConfig::EXP:
            os << "exponential decrease: " << search_config_.time_interval_ <<
                " times";
            break;
        case SLConfig::NO:
            os << "no limit";
            break;
    }

    os << "\nMaximum number of probes: " << search_config_.max_probes_number_;
    os << "\nProbes percentage: " << search_config_.probes_percentage_;
    os << "\nSplits: " << search_config_.splits_;
    os << "\nPenalty coefficient: " << search_config_.interval_penalty_;
    os << "\nLuby restarts scale: " << search_config_.luby_scale_;
    os << "\nRestart-on-fails probes: " << search_config_.fails_restart_probes_;
    os << "\nRestart-on-fails threshold: " << search_config_.fails_restart_thr_;
    os << "\nSynchronize with the validator: " << std::boolalpha <<
            search_config_.validator_synchronize_;
    os << "\nSubmit probes to the validator: " << std::boolalpha <<
            search_config_.submit_probes_;
    os << "\nRestart after interval: " << std::boolalpha <<
            search_config_.do_restarts_;
    os << '\n';

    return os;
}

void BalancingMonitor::ExitSearch() {
    initial_ss_size_ = 0;
}

void BalancingMonitor::BeginNextDecision(DecisionBuilder* const b) {
    if (initial_ss_size_ == 0) {
        initial_ss_size_ = CurrentSearchSpaceSize();
        paused_ = false;
    }

    if (!paused_ && sl_solver_.HelpAvailable()) {
        const auto ss_size = CurrentSearchSpaceSize();
        const double ss_part = double(ss_size) / initial_ss_size_;
        if (ss_part >= low_thr_ && ss_part <= high_thr_ && CanDetachSubTree()) {
            // Off-load some work
            // Take a snapshot of vars
            snapshot_asgn_.Store();

            // Log
            LOG4CXX_DEBUG(logger, "General balancer: Off-loading a region to a"
                    "helper: " << snapshot_asgn_.DebugString());

            CandidateVector work(1);
            FullAssignmentToLite(snapshot_asgn_, work.back().var_asgn_);
            sl_solver_.DispatchWork(work);

            // Fail will "detach" the tree from the current solver
            sl_solver_.BeginCustomFail();
            solver()->Fail();
        } else {
            // Cannot detach; try again later
            if (ss_part < low_thr_) {
                // Search space cannot increase until fail
                paused_ = true;
            }
        }
    }
}

void BalancingMonitor::BeginFail() {
    paused_ = false;
}

bool BalancingMonitor::CanDetachSubTree() const {
    /*
     *  First of all, we shouldn't be at the right-only path. Detaching
     *  at this point will end the current search.
     */
    if (solver()->SearchLeftDepth() == 0) {
        return false;
    }

    /*
     * Secondly, all domains must be without holes.
     */
    for (const IntVar *var: all_vars_) {
        if (var->Max() - var->Min() + 1 != var->Size()) {
            /*
             * We still might have a case of a variable x * c, which might
             * not have holes, but the above condition obviously would be
             * false.
             */
            auto hole_iter = var->MakeHoleIterator(true);
            hole_iter->Init();
            if (hole_iter->Ok()) {
                return false;
            }
        }
    }

    return true;
}

uint64_t BalancingMonitor::CurrentSearchSpaceSize() const {
    uint64_t ss_size = 1;
    for (const IntVar *var: all_vars_) {
        ss_size *= var->Size();
    }
    return ss_size;
}


} /* namespace searchlight */

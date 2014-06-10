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
 * This class represents a decision to set an interval for a variable.
 */
class SetIntervalDecision : public Decision {
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
    }

    /**
     * The right branch. It removes the interval.
     *
     * @param s the solver
     */
    virtual void Refute(Solver* const s) override {
        var_->RemoveInterval(min_, max_);
        // Temporarily deactivate the interval
        sl_search_.DeactivateVarInterval(var_, var_int_ind_, true);
        LOG4CXX_TRACE(logger, "Permanently deactivating interval for  " <<
                var_->DebugString() << ", [" << min_ << ", " << max_ << "]");
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
     * @param fails_limit the maximum number of fails before the search ends
     */
    IntervalImpactMonitor(Solver * const s, const int fails_limit) :
        SearchMonitor(s),
        fails_(0),
        leaves_(0),
        fails_limit_(fails_limit),
        problem_infeasible_(false) {}

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
     * Callback at the end of search
     */
    virtual void ExitSearch() override {
        if (solver()->state() == Solver::PROBLEM_INFEASIBLE) {
            problem_infeasible_ = true;
        }
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
    const int fails_limit_;

    // True, if the problem was infeasible (i.e., the interval is infeasible)
    bool problem_infeasible_;
};

/**
 * Search heuristic to randomly sample the specified interval for the
 * specified number of tries (fails).
 */
class IntervalImpactBuilder : public DecisionBuilder {
public:
    /**
     * Creates a new interval sample heuristic.
     *
     * @param search_vars decision variables
     * @param var the var to bind to an interval
     * @param min the left interval bound
     * @param max the right interval bound
     * @param monitors monitors to establish, vector will be modified inside,
     *     but left unchanged at exit
     */
    IntervalImpactBuilder(IntVarVector &search_vars, IntVar * const var,
            const int64_t min, const int64_t max,
            std::vector<SearchMonitor *> &monitors) :

            search_vars_(search_vars),
            var_(var), min_(min), max_(max),
            monitors_(monitors) {}

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

        // We are going to conduct a random search...
        DecisionBuilder * const random_db = s->MakePhase(search_vars_,
                Solver::CHOOSE_RANDOM, Solver::ASSIGN_RANDOM_VALUE);

        // with Luby restarts. TODO: change the scale factor?
        SearchMonitor * const luby_restart = s->MakeLubyRestart(1);

        // Nested search
        monitors_.push_back(luby_restart);
        s->Solve(random_db, monitors_);
        monitors_.pop_back();

        // After the nested search, we have nothing to do. Exit
        return nullptr;
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
    // All decision variables
    IntVarVector &search_vars_;

    // The variable to set
    IntVar * const var_;

    // The interval to set
    const int64_t min_, max_;

    // Monitors to establish for the explorer search
    std::vector<SearchMonitor *> &monitors_;
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
            }
            LOG4CXX_DEBUG(logger, "Finishing the search because "
                    "of a large fail ratio: " << ratio);

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

Decision* SLSearch::Next(Solver* const s) {
    if (!intervals_explored_) {
        InitIntervals(s, 1000);
        intervals_explored_ = true; // no rev on backtracking -- compute once
    }

    // Here, we have all impacts computed
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
        var_impacts_[var_ind].impacts_[int_ind].second.penalty_ /= 2;

        // and try to set it
        return s->RevAlloc(new SetIntervalDecision(*this, int_ind,
                var, int_min, int_max));
    }

    LOG4CXX_DEBUG(logger, "Initiating an in-interval search");

    // All variables are bound to intervals -- random search with restarts
    // We are going to conduct a random search...
    DecisionBuilder * const random_db = s->MakePhase(all_vars_,
            Solver::CHOOSE_RANDOM, Solver::ASSIGN_RANDOM_VALUE);

    // set the time limit if needed
    std::vector<SearchMonitor *> monitors(solver_montors_);
    if (search_time_limit_ >= 0) {
        LOG4CXX_DEBUG(logger, "In-interval search, time left: " <<
                search_time_limit_);
        monitors.push_back(s->MakeTimeLimit(search_time_limit_ / 2 * 1000));

        // Make Luby restarts for better coverage. TODO: play the scale factor?
        monitors.push_back(s->MakeLubyRestart(1));
    }
    /*
     * Establish a finish-on-frequent fails monitor. The rationale behind this
     * is that if we fail often before a leaf is reached, we are probably
     * wasting our time, since the search is random. If we don't have a time
     * limit, we might finish this search faster.
     *
     * Note, we don't establish Luby restarts if we don't have a time limit,
     * since the search will run infinitely long.
     */
    monitors.push_back(s->RevAlloc(new FinishOnFailsMonitor{s, 1000, 0.2}));

    // starting the timer
    const auto solve_start_time = std::chrono::steady_clock::now();
    s->Solve(random_db, monitors);
    // starting the timer
    const auto solve_end_time = std::chrono::steady_clock::now();

    // time elapsed (might be less than a half of search_time_limit_, if
    // FinishOnFails finished the search.
    int64_t solve_secodns =
            std::chrono::duration_cast<std::chrono::seconds>
            (solve_end_time - solve_start_time).count();
    LOG4CXX_DEBUG(logger, "The search went for " << solve_secodns << "s.");

    // Finished this combination of intervals, can safely fail
    if (search_time_limit_ > 0) {
        search_time_limit_ -= solve_secodns;
        if (search_time_limit_ < 0) {
            search_time_limit_ = 0;
        }
    }

    if (search_time_limit_ == 0) {
        LOG4CXX_DEBUG(logger, "In-interval search, no time left, terminating");
        dummy_monitor_.FinishCurrentSearch();
    } else {
        /*
         *  Will try a different interval combination. Since we penalize
         *  intervals we choose we should end up with a different combination,
         *  unless the current intervals are _reall_ promising (having a
         *  high impact).
         */
        LOG4CXX_DEBUG(logger, "In-interval search, restarting.");
        dummy_monitor_.RestartCurrentSearch();
    }

    return s->MakeFailDecision();
}

void SLSearch::InitIntervals(Solver * const s, const int64_t steps_limit) {
    LOG4CXX_DEBUG(logger, "Exploring interval impacts");
    std::vector<SearchMonitor *> nested_monitors(solver_montors_);
    for (auto &var_impact: var_impacts_) {
        for (auto &interval_impact: var_impact.impacts_) {
            const int64 min = interval_impact.first;
            const int64 max = min + var_impact.interval_length_ - 1;

            LOG4CXX_TRACE(logger, "Exploring interval [" <<
                    min << ", " << max << "] for " <<
                    var_impact.var_->DebugString());

            // TODO: play with the limit?
            IntervalImpactMonitor explorer_monitor{s, 1000};
            nested_monitors.push_back(&explorer_monitor);
            IntervalImpactBuilder explorer{all_vars_, var_impact.var_,
                    min, max, nested_monitors};
            s->Solve(&explorer);
            nested_monitors.pop_back();

            Impact &imp = interval_impact.second;
            imp.tried_assigns_ = explorer_monitor.GetTotalFails();
            imp.succ_assigns_ = explorer_monitor.GetLeavesReached();
            if (explorer_monitor.ProblemInfeasible()) {
                // We should deactivate the interval
                LOG4CXX_TRACE(logger, "The interval is infeasible");
                imp.active = false;
            }
            LOG4CXX_TRACE(logger, "Determined impact for " <<
                    "[" << min << ", " << max << "]: " << imp);
        }
    }
    LOG4CXX_DEBUG(logger, "Finished exploring interval impacts");
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


} /* namespace searchlight */

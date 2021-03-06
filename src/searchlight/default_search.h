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
 * @file default_search.h
 *
 * This is the default search heuristic used by Searchlight and the number
 * of utility classes that go with it.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_DEFAULT_SEARCH_H_
#define SEARCHLIGHT_DEFAULT_SEARCH_H_

#include "ortools_inc.h"
#include "base.h"
#include "searchlight.h"

namespace searchlight {

/**
 * This class represents the default searchlight search. This search basically
 * tries to first assign primary variables promising intervals and then explore
 * the intervals with the random search. Restarting, to cover more intervals,
 * is supported. The intervals are explored beforehand to assess their impact
 * by probing.
 */
class SLSearch : public DecisionBuilder {
public:
    /**
     * Constructs a new instance if the default searchlight search.
     *
     * @param sl the searchlight instance
     * @param sl_solver searchlight solver
     * @param primary_vars primary decision variables used for intervals
     * @param secondary_vars secondary decision variables used in the search
     * @param search_time_limit time limit in seconds; negative -- no limit
     */
    SLSearch(Searchlight &sl,
            SearchlightSolver &sl_solver,
            const IntVarVector &primary_vars,
            const IntVarVector &secondary_vars);

    /**
     * Produces the next decision.
     *
     * @param s the solver
     * @return next decision to branch or nullptr if at the leaf
     */
    virtual Decision* Next(Solver *s) override;

    /**
     * Returns a string for debugging.
     *
     * @return debug string
     */
    virtual std::string DebugString() const override {
        return "Searchlight Default DecisionBuilder";
    }

    /**
     * Outputs variable impacts to a stream.
     *
     * @param os stream for output
     * @return the same stream
     */
    std::ostream &OutputImpacts(std::ostream &os) const;

    /**
     * Outputs the search configuration.
     *
     * @param os stream to output to
     * @return the same stream
     */
    std::ostream &OutputConfig(std::ostream &os) const;

private:
    // Internal monitors and builders
    class IntervalImpactBuilder;
    class SetIntervalDecision;

    /*
     * Removes a variable's interval from the search process. After this it
     * will not be participating in the interval choosing process. The variable
     * must be registered at this decision builder.
     *
     * @param var the variable
     * @param int_num the ordinal number of the interval
     * @param reverse true, if the decision must be reversed at
     *     the backtracking restart; false if to remove permanently
     */
    void DeactivateVarInterval(const IntVar *var, int int_num, bool reverse) {
        assert(var_to_index_.find(var) != var_to_index_.end());

        const int var_ind = var_to_index_[var];
        bool &act_flag = var_impacts_[var_ind].impacts_[int_num].second.active;
        if (reverse) {
            solver_.SaveAndSetValue(&act_flag, false);
        } else {
            act_flag = false;
        }
    }

    // This is struct containing config params
    struct SLConfig {
        // How to spread time between intervals
        enum TimeStrategy {
            CONST, // constant
            EXP,   // exponential decrease
            NO,    // run each zoom-in region until exhaustion
        };
        TimeStrategy time_strategy_;

        // What strategy to use for the in-interval search
        enum InIntervalSearch {
            RANDOM,
            SPLIT
        };
        InIntervalSearch int_search_;

        // The length of the time interval (or down coeff. for EXP)
        int64_t time_interval_;

        // Maximum number of probes
        int max_probes_number_;

        // Percentage of probes for an interval
        double probes_percentage_;

        // Number of splits for each primarty vartiable's domain
        int splits_;

        // Penalty for the interval after it is chosen again
        double interval_penalty_;

        // Scale factor for Luby restarts
        int luby_scale_;

        // Periodicity of restart-on-fails monitor (in fails)
        int64_t fails_restart_probes_;

        // Leaves-to-fails threshold of restart for the monitor above
        double fails_restart_thr_;

        // Syncronize with the validator?
        bool validator_synchronize_;

        // Submit probes to the validator?
        bool submit_probes_;

        // Perform restarts after each interval?
        bool do_restarts_;

        SLConfig(const SearchlightConfig &sl_config) {
            const std::string time_strategy =
                    sl_config.get("searchlight.sl.time_strategy",
                            std::string("fin"));
            int64_t def_time_interval;
            if (time_strategy == "exp") {
                time_strategy_ = EXP;
                def_time_interval = 2;
            } else if (time_strategy == "const") {
                time_strategy_ = CONST;
                def_time_interval = 600; // 10 minutes
            } else if (time_strategy == "fin") {
                time_strategy_ = NO;
                def_time_interval = -1; // doesn't matter
            } else {
                throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                        SCIDB_LE_ILLEGAL_OPERATION) << "unknown time strategy";
            }

            const std::string int_search = sl_config.get(
                    "searchlight.sl.int_search",
                    std::string("split"));
            if (int_search == "random") {
                int_search_ = RANDOM;
            } else if (int_search == "split") {
                int_search_ = SPLIT;
            } else {
                throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                        SCIDB_LE_ILLEGAL_OPERATION)
                                << "unknown in-interval search";
            }

            time_interval_ = sl_config.get("searchlight.sl.time_interval",
                    def_time_interval);
            max_probes_number_ =
                    sl_config.get("searchlight.sl.max_probes", 1000);
            probes_percentage_ =
                    sl_config.get("searchlight.sl.probes_percentage", 0.01);
            splits_ =
                    sl_config.get("searchlight.sl.splits", 100);
            interval_penalty_ =
                    sl_config.get("searchlight.sl.interval_penalty", 0.5);

            luby_scale_ = sl_config.get("searchlight.sl.luby_scale", 1);

            fails_restart_probes_ =
                    sl_config.get("searchlight.sl.fails_restart_probes",
                            int64_t(1000));
            fails_restart_thr_ =
                    sl_config.get("searchlight.sl.fails_restart_thr", 0.2);
            validator_synchronize_ =
                    sl_config.get("searchlight.sl.val_sync", 1);

            submit_probes_ =
                    sl_config.get("searchlight.sl.probe_submit", 0);

            do_restarts_ =
                    sl_config.get("searchlight.sl.restarts", 0);
        }
    };
    // This describes information discovered for a variable's interval
    struct Impact {
        // Assignments tried for this interval
        int tried_assigns_ = 0;

        // Successful assignments discovered
        int succ_assigns_ = 0;

        // Approximate number of chunks would be read for the assignments
        int chunks_to_read_ = 0;

        // Penalty for the impact (used during dynamic search)
        double penalty_ = 1;

        // Is this impact used to assess the variable?
        bool active = true;

        // Returns "promising" impact: the number of successful assignments
        double PromiseImpact() const {
            double imp = tried_assigns_ == 0 ? 0.01 : // give it a chance
                    double(succ_assigns_) / tried_assigns_;
            return imp * penalty_;
        }
    };

    friend std::ostream &operator<<(std::ostream &s, const Impact &impact);

    // Interval impact: left point + the impact
    typedef std::pair<int64, Impact> IntervalImpact;

    // Describes a variable's interval impact
    struct VarImpactInfo {
        // The variable index
        int var_index = -1;

        // The ordinal of the interval
        int int_index = -1;

        // The impact
        double impact_ = -1;
    };

    struct IntervalVarImpacts {
        // The length of the interval for the variable (might be 1)
        const size_t interval_length_;

        // Variable impacts
        std::vector<IntervalImpact> impacts_;

        // The variable itself
        IntVar * const var_;

        IntervalVarImpacts(const size_t lbound, const size_t rbound,
                const size_t len, IntVar * const var) :
            interval_length_(len),
            var_(var) {
            // Empty impacts by default
            for (size_t i = lbound; i <= rbound; i += len) {
                impacts_.push_back(std::make_pair(i, Impact()));
            }
        }
    };

    /*
     * Explores all intervals. It will try steps_limit steps. A "step"
     * basically runs until a fail (leaf one or in-tree).
     */
    void InitIntervals(Solver *s);

    /*
     * Finds the best variable from the ones that are still interval unbound.
     * If all variables are bound to intervals, nullptr is returned. If a
     * valid variable is returned, min/max contain the best impact interval.
     */
    VarImpactInfo FindBestVar() const;

    /*
     * Checks if all primary variables are bound to intervals.
     */
    bool VarsIntervalBound() const;

    // Primary variables guiding the interval search (var --> index)
    std::unordered_map<const IntVar *, int> var_to_index_;

    // Impacts for the primary vars (indexes are in the map above)
    std::vector<IntervalVarImpacts> var_impacts_;

    // primary decision variables for impacts estimation
    const IntVarVector primary_vars_;

    // Decision variables not participating in discovering the intervals
    const IntVarVector secondary_vars_;

    // All decision variables combined for convenience
    IntVarVector all_vars_;

    // True, if we already explored the intervals
    bool intervals_explored_;

    // The searchlight instance
    Searchlight &sl_;

    // Searchlight solver
    SearchlightSolver &sl_solver_;

    // The main solver
    Solver &solver_;

    // SearchMonitor used to restart and finish the current search (hackish)
    SearchMonitor dummy_monitor_;

    // Search parameters
    const SLConfig search_config_;
};

/**
 * Balancing monitor implements solver balancing in the general case, when
 * there is no access to the search heuristic (decision builder).
 *
 * It monitors the state of the search space, and when it falls in the
 * specified interval, it detaches the corresponding sub-tree and sends it
 * to the helper. It then fails the search, thus cutting the sub-tree from
 * this solver.
 *
 * There is a number of restrictions for the balancer to work:
 *  1) It doesn't detach when the current path is all-right. The reason is that
 *      detaching at this point would end the current search, not just balance
 *      it.
 *  2) For a sub-tree to be detachable, all variables must not have holes
 *      in their domains. The reason behind that is to make serialization of
 *      domains faster. No holes means only intervals are required.
 *
 * User specifies threshold interval for possible off-loads. If the current
 * sub-tree is inside the interval, it can be detached. Otherwise,
 * the balancer waits until another sub-tree is available.
 */
class BalancingMonitor : public SearchMonitor {
public:
    /**
     * Creates a new balancing monitor.
     *
     * @param sl_solver searchlight solver to monitor
     * @param sl searchlight instance
     * @param all_vars all search variables
     * @param low_thr low threshold for candidate sub-tree
     * @param high_thr high threshold for candidate sub-tree
     */
    BalancingMonitor(SearchlightSolver &sl_solver, Searchlight &sl,
            const std::vector<IntVar *> &all_vars,
            double low_thr, double high_thr) :
                SearchMonitor{&sl_solver.GetSearchSolver()},
                low_thr_{low_thr},
                high_thr_{high_thr},
                all_vars_{all_vars},
                sl_(sl),
                sl_solver_(sl_solver),
                snapshot_asgn_(&sl_solver.GetSearchSolver()) {

        snapshot_asgn_.Add(all_vars_);
    }

    /**
     * Callback called at the end of the search.
     */
    virtual void ExitSearch() override;

    /**
     * Callback called before the decision builder produces a new decision.
     *
     * @param b decision builder
     */
    virtual void BeginNextDecision(DecisionBuilder* const b) override;

    /**
     * Callback called when the search starts failing, but before actual jump.
     */
    virtual void BeginFail() override;

private:
    // Compute current search space size
    uint64_t CurrentSearchSpaceSize() const;

    // Check if it is possible to detach current sub-tree
    bool CanDetachSubTree() const;

    // Low tree detach threshold
    const double low_thr_;

    // High tree detach threshold
    const double high_thr_;

    // Intial search space size
    uint64_t initial_ss_size_ = 0;

    // True, if the balancing is paused for now
    bool paused_ = false;

    // All variables in a single vector
    const std::vector<IntVar *> all_vars_;

    // Searchlight instance
    Searchlight &sl_;

    // SearchlightSolver
    SearchlightSolver &sl_solver_;

    // Assignment for taking snapshots
    Assignment snapshot_asgn_;
};


} /* namespace searchlight */
#endif /* SEARCHLIGHT_DEFAULT_SEARCH_H_ */

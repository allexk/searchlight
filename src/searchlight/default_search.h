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
     * @param solver the main solver
     * @param primary_vars primary decision variables used for intervals
     * @param secondary_vars secondary decision variables used in the search
     * @param splits the maximum number of intervals for a variable
     * @param search_time_limit time limit in seconds; negative -- no limit
     */
    SLSearch(const Searchlight &sl,
            Solver &solver,
            const IntVarVector &primary_vars,
            const IntVarVector &secondary_vars, size_t splits,
            int64_t search_time_limit) :
        secondary_vars_(secondary_vars),
        intervals_explored_(false),
        search_time_limit_(search_time_limit),
        sl_(sl),
        solver_(solver),
        solver_montors_(sl.GetMainMonitors()),
        dummy_monitor_(&solver) {

        for (auto int_var: primary_vars) {
            if (!int_var->Bound()) {
                // TODO: Do we really need this map?
                const int var_ind = var_to_index_.size();
                var_to_index_[int_var] = var_ind;

                const int64 var_min = int_var->Min();
                const int64 var_max = int_var->Max();
                const int64 var_len = var_max - var_min + 1;
                const size_t interval_len = var_len <= splits ?
                        1 : var_len / splits;
                var_impacts_.push_back(IntervalVarImpacts(var_min, var_max,
                        interval_len, int_var));
            }
        }

        all_vars_.reserve(primary_vars.size() + secondary_vars.size());
        all_vars_.insert(all_vars_.end(), primary_vars.begin(),
                primary_vars.end());
        all_vars_.insert(all_vars_.end(), secondary_vars.begin(),
                secondary_vars.end());
    }

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

private:
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
    void InitIntervals(Solver *s, int64_t steps_limit);

    /*
     * Finds the best variable from the ones that are still interval unbound.
     * If all variables are bound to intervals, nullptr is returned. If a
     * valid variable is returned, min/max contain the best impact interval.
     */
    VarImpactInfo FindBestVar() const;

    // Primary variables guiding the interval search (var --> index)
    std::unordered_map<const IntVar *, int> var_to_index_;

    // Impacts for the primary vars (indexes are in the map above)
    std::vector<IntervalVarImpacts> var_impacts_;

    // Decision variables not participating in discovering the intervals
    const IntVarVector secondary_vars_;

    // All decision variables combined for convenience
    IntVarVector all_vars_;

    // True, if we already explored the intervals
    bool intervals_explored_;

    // Time limit on the search in milliseconds (<0 -- no limit)
    int64_t search_time_limit_;

    // The searchlight instance
    const Searchlight &sl_;

    // The main solver
    Solver &solver_;

    // Monitors attached to the solver (at the searchlight)
    const std::vector<SearchMonitor *> &solver_montors_;

    // SearchMonitor used to restart and finish the current search (hackish)
    SearchMonitor dummy_monitor_;
};

} /* namespace searchlight */
#endif /* SEARCHLIGHT_DEFAULT_SEARCH_H_ */

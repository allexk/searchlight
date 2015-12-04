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
 * @file dist.cpp
 *
 * Implementation of the distance-based similarity search task.
 *
 * @author Alexander Kalinin
 */

#include <searchlight/searchlight.h>

#include <boost/lexical_cast.hpp>

namespace searchlight {

namespace {

/**
 * Creates an integer variable and assigns it a name.
 *
 * The deal here is that the solver may optimize the resulting variable by,
 * for example, altering its domain. It will still return a "correct" variable,
 * but the variable might be a cast variable for some internal expression.
 * Since names are really important in Searchlight for identification, we
 * enforce the name after the variable is created.
 *
 * @param solver solver to create the variable for
 * @param start left boundary for the domain (inclusive)
 * @param end right boundary for the domain (inclusive)
 * @param step step for values in the domain (1 means all values)
 * @param name the name for the variable
 * @return integer variable with the specified name
 */
IntVar *MakeIntVarWithName(Solver &solver, int32 start, int32 end,
        int32 step, const std::string &name) {
    IntVar *var;
    if (step == 1) {
        var = solver.MakeIntVar(start, end);
    } else {
        std::vector<int64> vals;
        for (int64 v = start; v <= end; v += step) {
            vals.push_back(v);
        }
        var = solver.MakeIntVar(vals);
    }
    var->set_name(name);
    return var;
}

}

extern "C"
void Dist(Searchlight *sl, uint32_t id) {
    // SL Solver
    SearchlightSolver &sl_solver = sl->GetSLSolver(id);
    // CP solver
    Solver &solver = sl_solver.GetSearchSolver();

    // or-tools is deterministic by default
    solver.ReSeed(ACMRandom::HostnamePidTimeSeed());

    // parse the params
    const SearchlightConfig &config = sl->GetConfig();

    // task params
    const int32 start_time = config.get<int32>("dist.l_time");
    const int32 end_time   = config.get<int32>("dist.u_time");
    const int32 query_dist = config.get<int32>("dist.dist", 0);
    const std::string query_seq = config.get("dist.query", "");
    const int32 step_time  = config.get<int32>("dist.step_time", 1);
    const int32 time_limit = config.get("dist.time_limit", 3600);
    const int luby_scale = config.get("searchlight.sl.luby_scale", 1);

    // problem variables
    std::vector<IntVar *> coords(1);
    coords[0] = MakeIntVarWithName(solver, start_time, end_time, step_time,
            "time");

    // attribute and sequence
    const std::string signal_name = config.get<std::string>("dist.signal");
    AttributeID attr = sl->RegisterAttribute(signal_name);
    AdapterPtr adapter = sl_solver.CreateAdapter("dist");
    const size_t seq_id = sl->RegisterQuerySequence(attr, query_seq);
    const size_t seq_len = sl->GetQuerySequence(seq_id).size();

    // valid sequence
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords[0], seq_len), end_time + 1));

    // create function
    UDFFunctionCreator sqdist_fab = sl->GetUDFFunctionCreator("sqdist");
    std::vector<int64> udf_params{attr, int64(seq_id), int64(seq_len), 0};
    IntExpr * const sqdist = solver.RevAlloc(sqdist_fab(&solver, adapter,
    		coords, udf_params));
    solver.AddConstraint(solver.MakeLessOrEqual(sqdist,
    		query_dist * query_dist));

    // create the search phase
    const std::string search_heuristic = config.get<std::string>("dist.db");
    DecisionBuilder *db;
    std::vector<SearchMonitor *> mons;
    if (search_heuristic == "impact") {
        db = solver.MakeDefaultPhase(coords);
    } else if (search_heuristic == "random") {
        db = solver.MakePhase(coords, Solver::CHOOSE_RANDOM,
            Solver::ASSIGN_RANDOM_VALUE);
        if (luby_scale != 0) {
            mons.push_back(solver.MakeLubyRestart(luby_scale));
        }
        mons.push_back(solver.MakeTimeLimit(time_limit * 1000));
    } else if (search_heuristic == "split") {
        db = solver.MakePhase(coords, Solver::CHOOSE_MAX_SIZE,
            Solver::SPLIT_LOWER_HALF);
        const bool balance = config.get("balance.solver_balance", 1);
        const double low_thr = config.get("balance.general_low", 0.1);
        const double high_thr = config.get("balance.general_high", 0.5);
        if (balance) {
            mons.push_back(sl_solver.CreateBalancingMonitor(coords,
                    low_thr, high_thr));
        }
    } else if (search_heuristic == "sl") {
        if (time_limit != 0) {
            mons.push_back(MakeCumulativeTimeLimit(solver, time_limit * 1000));
        }
        db = sl_solver.CreateDefaultHeuristic(coords, {});
    } else {
        db = solver.MakePhase(coords, Solver::CHOOSE_FIRST_UNBOUND,
            Solver::ASSIGN_MIN_VALUE);
    }

    // set task
    sl_solver.SetTask(coords, {}, db, mons);
}
} /* namespace searchlight */

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
 * @file sw.cpp
 *
 * Implementation of Semantic Windows (SW) task.
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
void MimicAvg(Searchlight *sl, uint32_t id) {
    // SL Solver
    SearchlightSolver &sl_solver = sl->GetSLSolver(id);
    // CP solver
    Solver &solver = sl_solver.GetSearchSolver();

    // or-tools is deterministic by default
    solver.ReSeed(ACMRandom::HostnamePidTimeSeed());

    // parse the params
    const SearchlightConfig &config = sl->GetConfig();

    // task params
    const int32 start_id = config.get<int32>("mimic.l_id");
    const int32 end_id   = config.get<int32>("mimic.u_id");
    const int32 start_time = config.get<int32>("mimic.l_time");
    const int32 end_time   = config.get<int32>("mimic.u_time");
    const int32 avg_l   = config.get<int32>("mimic.avg_l");
    const int32 avg_u   = config.get<int32>("mimic.avg_u");
    const int32 len_l  = config.get<int32>("mimic.len_l");
    const int32 len_u  = config.get<int32>("mimic.len_u");
    const int32 step_len  = config.get<int32>("mimic.step_len", 1);
    const int32 step_time  = config.get<int32>("mimic.step_time", 1);

    const int32 time_limit = config.get("mimic.time_limit", 3600);
    const int luby_scale = config.get("searchlight.sl.luby_scale", 1);

    // problem params
    std::vector<IntVar *> coords(2);
    std::vector<IntVar *> lens(2);

    // coords
    coords[0] = MakeIntVarWithName(solver, start_id, end_id, 1, "id");
    coords[1] = MakeIntVarWithName(solver, start_time, end_time, step_time,
            "time");

    // lens
    lens[0] = solver.MakeIntConst(1);
    lens[1] = MakeIntVarWithName(solver, len_l, len_u, step_len, "len");

    // convenience -- all vars in a single vector
    std::vector<IntVar *> all_vars(coords);
    all_vars.insert(all_vars.end(), lens.begin(), lens.end());

    // valid window
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords[1], lens[1]), end_time + 1));

    // average
    const std::string signal_name = config.get<std::string>("mimic.signal");
    AttributeID attr = sl->RegisterAttribute(signal_name);
    AdapterPtr adapter = sl_solver.CreateAdapter("mimic");
    UDFFunctionCreator avg_fab = sl->GetUDFFunctionCreator("avg");

    std::vector<int64> udf_params(1, int64(attr));
    IntExpr * const avg = solver.RevAlloc(avg_fab(&solver, adapter, all_vars,
            udf_params));
    solver.AddConstraint(solver.MakeBetweenCt(avg, avg_l, avg_u));

    // neighborhood
    const int32 lneighb_size  = config.get("mimic.neighborhood.l_size", 0);
    const int32 rneighb_size  = config.get("mimic.neighborhood.r_size", 0);
    if (lneighb_size || rneighb_size) {
        /*
         * FIXME: This is a mess! And should be rewritten ASAP!
         */
        UDFFunctionCreator max_fab = sl->GetUDFFunctionCreator("max");
        std::vector<IntVar *> part_vars(4);
        std::vector<IntVar *> part_maxs(2);

        // max of the region
        IntExpr * const region_max = solver.RevAlloc(max_fab(&solver, adapter,
                all_vars, udf_params));

        // left
        if (lneighb_size) {
            part_vars[0] = coords[0]; // the same id/timeline
            part_vars[1] = solver.MakeSum(coords[0], -lneighb_size)->Var();
            part_vars[2] = solver.MakeIntConst(1);
            part_vars[3] = solver.MakeIntConst(lneighb_size);

            // valid left neighborhood
            solver.AddConstraint(solver.MakeGreaterOrEqual(part_vars[1],
                    start_time));

            // difference of maximums
            const int32 left_max_diff  = config.get(
                    "mimic.neighborhood.left_max_diff", 0);
            IntVar * const left_max = solver.RevAlloc(max_fab(&solver,
                    adapter, part_vars, udf_params))->Var();
            solver.AddConstraint(solver.MakeGreater(solver.MakeDifference(
                    region_max, left_max), left_max_diff));
        }

        // right
        if (rneighb_size) {
            // right
            part_vars[0] = coords[0]; // the same id/timeline
            part_vars[1] = solver.MakeSum(coords[1], lens[1])->Var();
            part_vars[2] = solver.MakeIntConst(1);
            part_vars[3] = solver.MakeIntConst(rneighb_size);

            // valid right neighborhood
            solver.AddConstraint(solver.MakeLessOrEqual(
                    solver.MakeSum(part_vars[1], rneighb_size),
                    end_time + 1));

            // difference of maximums
            const int32 right_max_diff = config.get(
                    "mimic.neighborhood.right_max_diff", 0);
            IntVar * const right_max = solver.RevAlloc(max_fab(&solver,
                    adapter, part_vars, udf_params))->Var();
            solver.AddConstraint(solver.MakeGreater(solver.MakeDifference(
                    region_max, right_max), right_max_diff));
        }
    }

    // create the search phase
    const std::string search_heuristic = config.get<std::string>("mimic.db");
    DecisionBuilder *db;
    std::vector<SearchMonitor *> mons;
    if (search_heuristic == "impact") {
        db = solver.MakeDefaultPhase(all_vars);
    } else if (search_heuristic == "random") {
        db = solver.MakePhase(all_vars, Solver::CHOOSE_RANDOM,
            Solver::ASSIGN_RANDOM_VALUE);
        if (luby_scale != 0) {
            mons.push_back(solver.MakeLubyRestart(luby_scale));
        }
        mons.push_back(solver.MakeTimeLimit(time_limit * 1000));
    } else if (search_heuristic == "split") {
        db = solver.MakePhase(all_vars, Solver::CHOOSE_MAX_SIZE,
            Solver::SPLIT_LOWER_HALF);
        const bool balance = config.get("balance.solver_balance", 1);
        const double low_thr = config.get("balance.general_low", 0.1);
        const double high_thr = config.get("balance.general_high", 0.5);
        if (balance) {
            mons.push_back(sl_solver.CreateBalancingMonitor(all_vars, low_thr,
                    high_thr));
        }
    } else if (search_heuristic == "sl") {
        if (time_limit != 0) {
            mons.push_back(MakeCumulativeTimeLimit(solver, time_limit * 1000));
        }
        db = sl_solver.CreateDefaultHeuristic(coords, lens);
    } else {
        db = solver.MakePhase(all_vars, Solver::CHOOSE_FIRST_UNBOUND,
            Solver::ASSIGN_MIN_VALUE);
    }

    // set task
    sl_solver.SetTask(coords, lens, db, mons);
}
} /* namespace searchlight */

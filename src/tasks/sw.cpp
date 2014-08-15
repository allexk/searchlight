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
 * @file common.cpp
 *
 * Implementation of common taks function and classes.
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
void SemWindowsAvg(Searchlight *sl) {
    // get the solver
    Solver &solver = sl->GetSolver();

    // or-tools is deterministic by default
    solver.ReSeed(ACMRandom::HostnamePidTimeSeed());

    // parse the params
    const SearchlightConfig &config = sl->GetConfig();

    // task params
    const int32 start_x = config.get<int32>("sw.lx");
    const int32 end_x   = config.get<int32>("sw.ux");
    const int32 start_y = config.get<int32>("sw.ly");
    const int32 end_y   = config.get<int32>("sw.uy");
    const int32 avg_l   = config.get<int32>("sw.avg_l");
    const int32 avg_u   = config.get<int32>("sw.avg_u");
    const int32 size_l  = config.get<int32>("sw.size_l");
    const int32 size_u  = config.get<int32>("sw.size_u");
    const int32 len_lx  = config.get<int32>("sw.len_lx");
    const int32 len_ux  = config.get<int32>("sw.len_ux");
    const int32 len_ly  = config.get<int32>("sw.len_ly");
    const int32 len_uy  = config.get<int32>("sw.len_uy");
    const int32 step_x  = config.get<int32>("sw.step_x", 1);
    const int32 step_y  = config.get<int32>("sw.step_y", 1);
    const int32 step_lx  = config.get<int32>("sw.step_lx", 1);
    const int32 step_ly  = config.get<int32>("sw.step_ly", 1);

    const int32 time_limit = config.get("sw.time_limit", 3600);
    const int luby_scale = config.get("searchlight.sl.luby_scale", 1);
    const int splits = config.get("sw.splits", 100);

    // problem params
    std::vector<IntVar *> coords(2);
    std::vector<IntVar *> lens(2);

    // coords
    coords[0] = MakeIntVarWithName(solver, start_x, end_x, step_x, "x");
    coords[1] = MakeIntVarWithName(solver, start_y, end_y, step_y, "y");

    // lens
    lens[0] = MakeIntVarWithName(solver, len_lx, len_ux, step_lx, "lx");
    lens[1] = MakeIntVarWithName(solver, len_ly, len_uy, step_ly, "ly");

    // convenience -- all vars in a single vector
    std::vector<IntVar *> all_vars(coords);
    all_vars.insert(all_vars.end(), lens.begin(), lens.end());

    // valid window
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords[0], lens[0]), end_x + 1));
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords[1], lens[1]), end_y + 1));
    solver.AddConstraint(solver.MakeBetweenCt(
            solver.MakeProd(lens[0], lens[1]), size_l, size_u));

    // average
    AttributeID attr = sl->RegisterAttribute("val");
    AdapterPtr adapter = sl->CreateAdapter("sw");
    UDFFunctionCreator avg_fab = sl->GetUDFFunctionCreator("avg");

    std::vector<int64> udf_params(1, int64(attr));
    IntExpr * const avg = solver.RevAlloc(avg_fab(&solver, adapter, all_vars,
            udf_params));
    solver.AddConstraint(solver.MakeBetweenCt(avg, avg_l, avg_u));

    // neighborhood
    const int32 neighb_size  = config.get("sw.neighborhood.size", 0);
    if (neighb_size) {
        /*
         * FIXME: This is a mess! And should be rewritten ASAP!
         */
        UDFFunctionCreator max_fab = sl->GetUDFFunctionCreator("max");
        UDFFunctionCreator min_fab = sl->GetUDFFunctionCreator("min");
        std::vector<IntVar *> part_vars(4);
        std::vector<IntVar *> part_maxs(4);

        // left
        part_vars[0] = solver.MakeSum(coords[0], -neighb_size)->Var();
        part_vars[1] = coords[1];
        part_vars[2] = solver.MakeIntConst(neighb_size);
        part_vars[3] = lens[1];
        part_maxs[0] = solver.RevAlloc(max_fab(&solver,
                adapter, part_vars, udf_params))->Var();

        // top
        part_vars[0] = solver.MakeSum(coords[0], -neighb_size)->Var();
        part_vars[1] = solver.MakeSum(coords[1], -neighb_size)->Var();
        part_vars[2] = solver.MakeSum(lens[0], 2 * neighb_size)->Var();
        part_vars[3] = solver.MakeIntConst(neighb_size);
        part_maxs[1] = solver.RevAlloc(max_fab(&solver,
                adapter, part_vars, udf_params))->Var();

        // right
        part_vars[0] = solver.MakeSum(coords[0], lens[0])->Var();
        part_vars[1] = coords[1];
        part_vars[2] = solver.MakeIntConst(neighb_size);
        part_vars[3] = lens[1];
        part_maxs[2] = solver.RevAlloc(max_fab(&solver,
                adapter, part_vars, udf_params))->Var();

        // bottom
        part_vars[0] = solver.MakeSum(coords[0], -neighb_size)->Var();
        part_vars[1] = solver.MakeSum(coords[1], lens[1])->Var();
        part_vars[2] = solver.MakeSum(lens[0], 2 * neighb_size)->Var();
        part_vars[3] = solver.MakeIntConst(neighb_size);
        part_maxs[3] = solver.RevAlloc(max_fab(&solver,
                adapter, part_vars, udf_params))->Var();

        // max of the neighborhood
        IntExpr * const nmax = solver.MakeMax(part_maxs);

        const int32 max_diff  = config.get("sw.neighborhood.max_diff", 0);
        if (max_diff) {
            // max of the region
            IntExpr * const rmax = solver.RevAlloc(max_fab(&solver, adapter,
                    all_vars, udf_params));
            solver.AddConstraint(solver.MakeGreater(
                    solver.MakeDifference(rmax, nmax), max_diff));
        }

        const int32 min_diff  = config.get("sw.neighborhood.min_diff", 0);
        if (min_diff) {
            // min of the region
            IntExpr * const rmin = solver.RevAlloc(min_fab(&solver, adapter,
                    all_vars, udf_params));
            solver.AddConstraint(solver.MakeGreater(
                    solver.MakeDifference(rmin, nmax), min_diff));
        }
    }

    // create the search phase
    const std::string search_heuristic = config.get<std::string>("sw.db");
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
    } else if (search_heuristic == "sl") {
        if (time_limit != 0) {
            mons.push_back(MakeCumulativeTimeLimit(solver, time_limit * 1000));
        }
        db = solver.RevAlloc(
                sl->CreateDefaultHeuristic(coords, lens, splits));
    } else {
        db = solver.MakePhase(all_vars, Solver::CHOOSE_FIRST_UNBOUND,
            Solver::ASSIGN_MIN_VALUE);
    }

    // solve!
    sl->Prepare(coords, lens, db, mons);
}
} /* namespace searchlight */

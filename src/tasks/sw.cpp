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

extern "C"
void SemWindowsAvg(Searchlight *sl) {
    // get the solver
    Solver &solver = sl->GetSolver();

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

    const int32 time_limit = config.get("sw.time_limit", 3600);
    const int luby_scale = config.get("searchlight.sl", 1);
    const int splits = config.get("sw.splits", 100);

    // problem params
    std::vector<IntVar *> coords(2);
    std::vector<IntVar *> lens(2);
    coords[0] = solver.MakeIntVar(start_x, end_x, "x");
    coords[1] = solver.MakeIntVar(start_y, end_y, "y");
    lens[0] = solver.MakeIntVar(len_lx, len_ux, "lx");
    lens[1] = solver.MakeIntVar(len_ly, len_uy, "ly");

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
    IntExpr *avg = solver.RevAlloc(avg_fab(&solver, adapter, all_vars,
            udf_params));
    solver.AddConstraint(solver.MakeBetweenCt(avg, avg_l, avg_u));

    // create the search phase
    const std::string search_heuristic = config.get<std::string>("sw.db");
    DecisionBuilder *db;
    std::vector<SearchMonitor *> mons;
    if (search_heuristic == "impact") {
        db = solver.MakeDefaultPhase(all_vars);
    } else if (search_heuristic == "random") {
        db = solver.MakePhase(all_vars, Solver::CHOOSE_RANDOM,
            Solver::ASSIGN_RANDOM_VALUE);
        mons.push_back(solver.MakeLubyRestart(luby_scale));
        mons.push_back(solver.MakeTimeLimit(time_limit * 1000));
    } else if (search_heuristic == "split") {
        db = solver.MakePhase(all_vars, Solver::CHOOSE_MAX_SIZE,
            Solver::SPLIT_LOWER_HALF);
    } else if (search_heuristic == "sl") {
        db = solver.RevAlloc(
                sl->CreateDefaultHeuristic(coords, lens, splits, time_limit));
    } else {
        db = solver.MakePhase(all_vars, Solver::CHOOSE_FIRST_UNBOUND,
            Solver::ASSIGN_MIN_VALUE);
    }

    // solve!
    sl->Solve(db, all_vars, mons);
}
} /* namespace searchlight */

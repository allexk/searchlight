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

#include "common.h"
#include <searchlight/searchlight.h>

#include <boost/lexical_cast.hpp>

namespace searchlight {

extern "C"
void SemWindowsAvg(Searchlight *sl, const std::string &params) {
    // get the solver
    Solver &solver = sl->GetSolver();

    // parse the params
    ParamsMap pmap;
    ParseParameters(params, pmap);

    // task params
    const int32 start_x = boost::lexical_cast<int32>(pmap["lx"]);
    const int32 end_x = boost::lexical_cast<int32>(pmap["ux"]);
    const int32 start_y = boost::lexical_cast<int32>(pmap["ly"]);
    const int32 end_y = boost::lexical_cast<int32>(pmap["uy"]);
    const int32 avg_l = boost::lexical_cast<int32>(pmap["avg_l"]);
    const int32 avg_u = boost::lexical_cast<int32>(pmap["avg_u"]);
    const int32 size_l = boost::lexical_cast<int32>(pmap["size_l"]);
    const int32 size_u = boost::lexical_cast<int32>(pmap["size_u"]);

    // problem params
    std::vector<IntVar *> coords(2);
    std::vector<IntVar *> lens(2);
    coords[0] = solver.MakeIntVar(start_x, end_x, "x");
    coords[1] = solver.MakeIntVar(start_y, end_y, "y");
    lens[0] = solver.MakeIntVar(1, size_u, "lx");
    lens[1] = solver.MakeIntVar(1, size_u, "ly");

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
    DecisionBuilder * const db = solver.MakePhase(all_vars,
            Solver::CHOOSE_FIRST_UNBOUND,
            Solver::ASSIGN_MIN_VALUE);

    // solve!
    std::vector<SearchMonitor *> mons;
    sl->Solve(db, all_vars, mons);
}
} /* namespace searchlight */

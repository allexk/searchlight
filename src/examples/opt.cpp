/*
 * Copyright 2013, Brown University, Providence, RI.
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
 * @file sum.cpp
 * This is an example of using an aggregate CSP on a 2D array.
 *
 * @author Alexander Kalinin
 */
#include "gflags/gflags.h"
#include "constraint_solver/constraint_solveri.h"
#include "cpdb/aggr.h"
#include "cpdb/simple_array_2d.h"
#include <fstream>

DEFINE_string(array_file, "", "A file to import the array from.");
DEFINE_bool(print_array, false, "Dump the array before solving");
DEFINE_bool(mm_array, true, "True if the array is in MatrixMarket format");

namespace cpdb {

using operations_research::Solver;
using operations_research::IntVar;
using operations_research::IntExpr;
using operations_research::SolutionCollector;
using operations_research::DecisionBuilder;
using operations_research::OptimizeVar;

void SolveOptProblem(const SimpleArray2D &array) {
    Solver solver("array opt problem");

    // we treat columns as x-axis and rows as y-axis
    const int32 end_x = array.ColsNum() - 1;
    const int32 end_y = array.RowsNum() - 1;

    // problem params
    std::vector<IntVar *> coords(2);
    std::vector<IntVar *> lens(2);
    coords[0] = solver.MakeIntVar(0, end_x, "x");
    coords[1] = solver.MakeIntVar(0, end_y, "y");
    lens[0] = solver.MakeIntVar(3, 3, "lx");
    lens[1] = solver.MakeIntVar(2, 2, "ly");

    // valid window
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords[0], lens[0]), end_x + 1));
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords[1], lens[1]), end_y + 1));

    // convenience -- all vars in a single vector
    std::vector<IntVar *> all_vars(coords);
    all_vars.insert(all_vars.end(), lens.begin(), lens.end());

    // create the sum expression
    AggrFuncExpr * const sum_expr = solver.RevAlloc(new AggrFuncExpr(&solver,
            array, SUM, coords, lens));
    IntVar *sum_var = sum_expr->Var();
    sum_var = solver.RegisterIntVar(sum_var);

    // add the main constraint
    OptimizeVar * const opt = solver.MakeMaximize(sum_var, 1);

    // create the solutions collector
    SolutionCollector * const sols = solver.MakeLastSolutionCollector();
    sols->Add(all_vars);
    sols->Add(sum_var);

    // create the search phase
    DecisionBuilder * const db = solver.MakePhase(all_vars,
            Solver::CHOOSE_FIRST_UNBOUND,
            Solver::ASSIGN_MIN_VALUE);

    // solve!
    solver.Solve(db, opt, sols);

    // print solutions
    const int sol_num = sols->solution_count();
    for (int i = 0; i < sol_num; i++) {
        printf("Solution[%d]:\n", i);
        printf("\tx = %d,y = %d, lx = %d, ly = %d, sum = %d\n",
                sols->Value(i, coords[0]),
                sols->Value(i, coords[1]),
                sols->Value(i, lens[0]),
                sols->Value(i, lens[1]),
                sols->Value(i, sum_var));
    }

    // print stats
    sum_expr->PrintStats();
}
}

int main(int argc, char **argv) {
    google::ParseCommandLineFlags(&argc, &argv, true);

    cpdb::SimpleArray2D array;
    array.LoadFromFile(FLAGS_array_file.c_str(), FLAGS_mm_array);
    if (FLAGS_print_array) {
        std::ofstream dump_file("aray.dump");
        array.DumpArray(dump_file);
    }

    cpdb::SolveOptProblem(array);

    return 0;
}

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
 * @file sim.cpp
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

void SolveSimProblem(const SimpleArray2D &array) {
    Solver solver("array sim problem");

    // we treat columns as x-axis and rows as y-axis
    const int32 end_x = array.ColsNum() - 1;
    const int32 end_y = array.RowsNum() - 1;

    // problem params:
    // first window
    std::vector<IntVar *> coords1(2);
    std::vector<IntVar *> lens1(2);
    coords1[0] = solver.MakeIntVar(0, end_x, "x");
    coords1[1] = solver.MakeIntVar(0, end_y, "y");
    lens1[0] = solver.MakeIntVar(2, 2, "lx");
    lens1[1] = solver.MakeIntVar(2, 2, "ly");
    // validity
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords1[0], lens1[0]), end_x + 1));
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords1[1], lens1[1]), end_y + 1));

    // second window
    std::vector<IntVar *> coords2(2);
    std::vector<IntVar *> lens2(2);
    coords2[0] = solver.MakeIntVar(0, end_x, "x");
    coords2[1] = solver.MakeIntVar(0, end_y, "y");
    lens2[0] = solver.MakeIntVar(2, 2, "lx");
    lens2[1] = solver.MakeIntVar(2, 2, "ly");
    // validity
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords2[0], lens2[0]), end_x + 1));
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeSum(coords2[1], lens2[1]), end_y + 1));

    // convenience -- all vars in a single vector
    std::vector<IntVar *> all_vars(coords1);
    all_vars.insert(all_vars.end(), lens1.begin(), lens1.end());
    all_vars.insert(all_vars.end(), coords2.begin(), coords2.end());
    all_vars.insert(all_vars.end(), lens2.begin(), lens2.end());

    // create the sum expression
    // first window
    AggrFuncExpr * const sum_expr1 = solver.RevAlloc(new AggrFuncExpr(&solver,
            array, SUM, coords1, lens1));
    IntVar *sum_var1 = sum_expr1->Var();
    sum_var1 = solver.RegisterIntVar(sum_var1);

    // second window
    AggrFuncExpr * const sum_expr2 = solver.RevAlloc(new AggrFuncExpr(&solver,
            array, SUM, coords2, lens2));
    IntVar *sum_var2 = sum_expr2->Var();
    sum_var2 = solver.RegisterIntVar(sum_var2);

    // add constraint on the order --- to break the symmetry
    // rank is a major-row number of the corner
    IntExpr *rank1 = solver.MakeSum(solver.MakeProd(coords1[1], end_x + 1),
            coords1[0]);
    IntExpr *rank2 = solver.MakeSum(solver.MakeProd(coords2[1], end_x + 1),
            coords2[0]);
    solver.AddConstraint(solver.MakeLess(rank1, rank2));

    // add the main constraints: on a sum and between sums
    //sump.AddConstraint(sump.MakeLess(sum_var, 10));
    solver.AddConstraint(solver.MakeBetweenCt(sum_var1, 5, 10));
    //solver.AddConstraint(solver.MakeEquality(sum_expr1, sum_expr2));
    solver.AddConstraint(solver.MakeLessOrEqual(
            solver.MakeAbs(solver.MakeDifference(sum_expr1, sum_expr2)), 1));

    // create the solutions collector
    SolutionCollector * const sols = solver.MakeAllSolutionCollector();
    sols->Add(all_vars);
    sols->Add(sum_var1);
    sols->Add(sum_var2);

    // create the search phase
    DecisionBuilder * const db = solver.MakePhase(all_vars,
            Solver::CHOOSE_FIRST_UNBOUND,
            Solver::ASSIGN_MIN_VALUE);

    // solve!
    solver.Solve(db, sols);

    // print solutions
    const int sol_num = sols->solution_count();
    for (int i = 0; i < sol_num; i++) {
        printf("Solution[%d]:\n", i);
        printf("\tFirst: x = %d,y = %d, lx = %d, ly = %d, sum = %d\n",
                sols->Value(i, coords1[0]),
                sols->Value(i, coords1[1]),
                sols->Value(i, lens1[0]),
                sols->Value(i, lens1[1]),
                sols->Value(i, sum_var1));
        printf("\tSecond: x = %d,y = %d, lx = %d, ly = %d, sum = %d\n",
                sols->Value(i, coords2[0]),
                sols->Value(i, coords2[1]),
                sols->Value(i, lens2[0]),
                sols->Value(i, lens2[1]),
                sols->Value(i, sum_var2));
    }

    // print stats
    sum_expr1->PrintStats();
    sum_expr2->PrintStats();
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

    cpdb::SolveSimProblem(array);

    return 0;
}

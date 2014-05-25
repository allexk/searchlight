/*
 * Copyright 2014, Brown University, Providence, RI.
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
 * @file searchlight_collector.cpp
 * The implementation of the searchlight solution collector.
 *
 * @author Alexander Kalinin
 */

#include "searchlight_collector.h"
#include "searchlight_task.h"

namespace searchlight {

TaskSolutionCollector::TaskSolutionCollector(SearchlightTask &task, Solver *s) :
    SolutionCollector(s),
    task_(task) {}

void TaskSolutionCollector::ExitSearch() {
    task_.OnFinishSearch();
}

bool TaskSolutionCollector::AtSolution() {
    // we can reuse the same prototype and do not have to store the values
    const Assignment::IntContainer &vars = prototype_.get()->IntVarContainer();

    // stringify the solution
    std::ostringstream sol_string;
    for (size_t i = 0; i < vars.Size(); i++) {
        const IntVar *v = vars.Element(i).Var();
        sol_string << v->name() << "=" << v->Value();
        if (i != vars.Size() - 1) {
            sol_string << ", ";
        }
    }

    task_.AddSolution(sol_string.str());
    return true;
}
}

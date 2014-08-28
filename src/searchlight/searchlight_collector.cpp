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

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.collector"));

// The logger for results
static log4cxx::LoggerPtr result_logger(
        log4cxx::Logger::getLogger("searchlight.result"));

TaskSolutionCollector::TaskSolutionCollector(SearchlightTask &task, Solver *s) :
    SolutionCollector(s),
    task_(task) {}

bool TaskSolutionCollector::AtSolution() {
    // we can reuse the same prototype and do not have to store the values
    const Assignment::IntContainer &vars = prototype_.get()->IntVarContainer();

    // stringify the solution
    std::vector<int64_t> vals(vars.Size());
    for (size_t i = 0; i < vars.Size(); i++) {
        const IntVar *v = vars.Element(i).Var();
        vals[i] = v->Value();
    }
    task_.ReportSolution(vals);

    return true;
}

std::string TaskSolutionCollector::SolutionToString(
        const std::vector<int64_t> &vals) const {
    // we can reuse the same prototype and do not have to store the values
    const Assignment::IntContainer &vars = prototype_.get()->IntVarContainer();

    // stringify the solution
    std::ostringstream sol_string;
    for (size_t i = 0; i < vars.Size(); i++) {
        const IntVar *v = vars.Element(i).Var();
        sol_string << v->name() << "=" << vals[i];
        if (i != vars.Size() - 1) {
            sol_string << ", ";
        }
    }
    const std::string solution = sol_string.str();
    LOG4CXX_TRACE(logger, "Collected a new solution: " << solution);
    LOG4CXX_TRACE(result_logger, solution);

    return solution;
}

}

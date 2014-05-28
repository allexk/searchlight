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
 * @file searchlight.cpp
 * The implementation of the main searchlight module.
 *
 * @author Alexander Kalinin
 */

#include "searchlight.h"
#include "validator.h"
#include "searchlight_collector.h"

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.searchlight"));

Searchlight::~Searchlight() {
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(
            total_solve_time_).count();
    const auto usecs = std::chrono::duration_cast<std::chrono::milliseconds>(
            total_solve_time_).count();
    LOG4CXX_INFO(logger, "Searchlight solve time: " << secs << '.' <<
            usecs << 's');

    delete array_desc_;
}

bool Searchlight::Solve(DecisionBuilder *db, const IntVarVector &vars,
        const std::vector<SearchMonitor *> &monitors) {
    /*
     * First, we need to establish our own validation collector, create
     * a validator and pass it the names of the decision vars
     */
    StringVector var_names(vars.size());
    int i = 0;
    for (IntVarVector::const_iterator cit = vars.begin(); cit != vars.end();
            cit++) {
        const IntVar * const var = *cit;
        if (!var->HasName()) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Every decision variable must have a name!";
        }
        var_names[i++] = var->name();
    }

    // establish the validator
    if (!collector_) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "No solution collector registered!";
    }
    LOG4CXX_INFO(logger, "Initiating the validator");
    Validator validator(*this, var_names, *collector_);
    boost::thread validator_thread(boost::ref(validator));

    // Establish monitors: validator (to transfer leaves) and terminator
    std::vector<SearchMonitor *> solver_monitors(monitors);
    ValidatorMonitor val_monitor(validator, vars, &solver_);
    SearchLimit *terminator = solver_.MakeCustomLimit(
            NewPermanentCallback(this, &Searchlight::CheckTerminate));
    solver_monitors.push_back(&val_monitor);
    solver_monitors.push_back(terminator);

    // starting the timer
    const auto solve_start_time = std::chrono::steady_clock::now();

    // start the search
    LOG4CXX_INFO(logger, "Starting the main search");
    solver_.Solve(db, solver_monitors);

    // Terminate validator
    LOG4CXX_INFO(logger, "Signaling the validator and waiting");
    validator.SignalEnd();
    validator_thread.join();

    // stopping the timer
    const auto solve_end_time = std::chrono::steady_clock::now();
    total_solve_time_ =
            std::chrono::duration_cast<decltype(total_solve_time_)>(
            solve_end_time - solve_start_time);

    LOG4CXX_INFO(logger, "Finished the main search");
    return validator.GetValidatorSolverResult();
}

bool ValidatorMonitor::AtSolution() {
    Assignment * const asgn = prototype_.get();

    // Store the solution (assuming complete assignment)
    asgn->Store();
    LOG4CXX_DEBUG(logger, "Encountered a leaf: " << asgn->DebugString());

    // Should check if we have a complete assignment
    bool complete = true;
    for (IntVarVector::const_iterator cit = vars_.begin(); cit != vars_.end();
            cit++) {
        const IntVar * const var = *cit;
        if (!asgn->Bound(var)) {
            complete = false;
            break;
        }
    }

    if (complete) {
        validator_.AddSolution(*asgn);
    }

    return true;
}
} /* namespace searchlight */

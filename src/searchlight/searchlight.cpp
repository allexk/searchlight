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
#include "default_search.h"

/*
 * FIXME:
 * This is a local patched copy of the corresponding Boost file. It fixes a
 * well-known bug on C++11. Since I do not have a newer version of Boost on
 * the development machine, thus the hack. The bug is fixed upstream.
 */
//#include "json_parser_read.hpp"

/*
 * FIXME: There is a compile warning here about an always-true comparison in
 * json_parser_write inside the header. Actually, the warning is
 * valid for char-based strings, but harmless. See boost bug #5598
 */
#include <boost/property_tree/json_parser.hpp>

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.searchlight"));

Searchlight::~Searchlight() {
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(
            total_solve_time_).count();
    const auto usecs = std::chrono::duration_cast<std::chrono::milliseconds>(
            total_solve_time_).count();
    LOG4CXX_INFO(logger, "Solver total time: " << secs << '.' <<
            usecs << 's');

    search_monitors_.Disband();
    EndAndDestroyValidator();

    delete array_desc_;
}

void Searchlight::Prepare(const IntVarVector &primary_vars,
        const IntVarVector &secondary_vars,
        DecisionBuilder *db, const std::vector<SearchMonitor *> &monitors) {
    /*
     * Perform initialization
     */
    primary_vars_ = primary_vars;
    secondary_vars_ = secondary_vars;
    db_ = db;

    /*
     * First, we need to establish our own validation collector, create
     * a validator and pass it the names of the decision vars
     */
    IntVarVector vars{primary_vars};
    vars.insert(vars.end(), secondary_vars.begin(), secondary_vars.end());

    StringVector var_names(vars.size());
    int i = 0;
    for (IntVarVector::const_iterator cit = vars.begin(); cit != vars.end();
            cit++) {
        const IntVar * const var = *cit;
        std::string var_name;
        if (var->HasName()) {
            var_name = var->name();
        } else {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                    SCIDB_LE_ILLEGAL_OPERATION)
                    << "Every decision variable must have a name!";
        }
        var_names[i++] = var_name;
    }

    /*
     * Determine our partition. First, determine the variable with the largest
     * domain and the assume equivalent stripes along its domain.
     */
    IntVar *split_var = nullptr;
    uint64_t max_size = 0;
    for (const auto var: primary_vars) {
        const uint64_t var_size = var->Max() - var->Min() + 1;
        if (!split_var || var_size > max_size) {
            split_var = var;
            max_size = var_size;
        }
    }
    assert(split_var);

    // define partition
    const int64_t partition_len = max_size / instance_count_;
    const int64_t partition_start = split_var->Min() +
            partition_len * instance_id_;
    const int64_t partition_end = partition_start + partition_len - 1;
    split_var->SetRange(partition_start, partition_end);
    LOG4CXX_INFO(logger, "Set partition for var=" << split_var->DebugString());

    // logging
    if (logger->isDebugEnabled()) {
        /*
         *  Unfortunately, it will dump it into std::cerr, which is
         *  redirected by SciDb to a file.
         */
        ModelVisitor *pmv = solver_.MakePrintModelVisitor();
        solver_.Accept(pmv);
    }

    // establish the validator
    if (!search_monitors_.collector_) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "No solution collector registered!";
    }
    LOG4CXX_INFO(logger, "Initiating the validator");
    validator_ = new Validator(*this, var_names, *search_monitors_.collector_);

    // Establish monitors: validator (to transfer leaves) and terminator
    search_monitors_.user_monitors_ = monitors;
    search_monitors_.validator_monitor_ =
            new ValidatorMonitor(*validator_, vars, &solver_);
    SearchLimit *terminator = solver_.MakeCustomLimit(
            NewPermanentCallback(this, &Searchlight::CheckTerminate));
    search_monitors_.aux_monitors_.push_back(terminator);
}

void Searchlight::Solve() {
    if (!validator_thread_) {
        LOG4CXX_INFO(logger, "Starting the validator thread");
        validator_thread_ = new std::thread(std::ref(*validator_));
    }

    // Starting the timer
    const auto solve_start_time = std::chrono::steady_clock::now();

    // start the search
    LOG4CXX_INFO(logger, "Starting the main search");
    std::vector<SearchMonitor *> solve_monitors(
            search_monitors_.user_monitors_);
    solve_monitors.insert(solve_monitors.end(),
            search_monitors_.aux_monitors_.begin(),
            search_monitors_.aux_monitors_.end());
    solve_monitors.push_back(search_monitors_.validator_monitor_);
    solver_.Solve(db_, solve_monitors);

    // stopping the timer
    const auto solve_end_time = std::chrono::steady_clock::now();
    total_solve_time_ =
            std::chrono::duration_cast<decltype(total_solve_time_)>(
            solve_end_time - solve_start_time);

    /*
     * Output some stats. Note that the solver's fail counter contains
     * not only true fails (due to violations), but also in-solution fails,
     * when we it just backtracks to find the next one. In our case, we always
     * continue after another candidate is found, so fails-sols should give
     * the number of true fails.
     */
    int64 total_fails = solver_.failures();
    int64_t total_candidates =
            search_monitors_.validator_monitor_->CandidatesNumber();
    LOG4CXX_INFO(logger, "Main search stats: fails=" << total_fails <<
            ", true fails=" << (total_fails - total_candidates) <<
            ", candidates=" << total_candidates);

    LOG4CXX_INFO(logger, "Finished the main search");
    // TODO: temporary; need a better solution when balancing comes into play
    validator_->SignalEnd();
}

void Searchlight::EndAndDestroyValidator() {
    if (validator_) {
        if (validator_thread_) {
            LOG4CXX_INFO(logger, "Signaling the validator and waiting");
            validator_->SignalEnd();
            validator_thread_->join();
            delete validator_thread_;
        }

        delete validator_;
    }
}

bool ValidatorMonitor::AtSolution() {
    Assignment * const asgn = prototype_.get();

    // Store the solution (assuming complete assignment)
    asgn->Store();
    LOG4CXX_TRACE(logger, "Encountered a leaf: " << asgn->DebugString() <<
            ", depth=" << solver()->SearchDepth());

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
        candidates_++;
        validator_.AddSolution(*asgn);
    }

    return true;
}

DecisionBuilder *Searchlight::CreateDefaultHeuristic(
        const IntVarVector &primary_vars,
        const IntVarVector &secondary_vars, size_t splits) {
    return new SLSearch(*this, solver_, primary_vars, secondary_vars, splits);
}

void Searchlight::ReadConfig(const std::string &file_name) {
    try {
        boost::property_tree::read_json(file_name, config_);
    } catch (const boost::property_tree::json_parser_error &e) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << e.what();
    }

    if (logger->isDebugEnabled()) {
        std::ostringstream deb;
        deb << "Config for the task follows:\n";
        boost::property_tree::write_json(deb, config_, true);
        logger->debug(deb.str(), LOG4CXX_LOCATION);
    }
}

SearchMonitor *MakeCumulativeTimeLimit(Solver &s, int64 time_ms) {
    return s.MakeLimit(time_ms, kint64max, kint64max, kint64max, true, true);
}

} /* namespace searchlight */

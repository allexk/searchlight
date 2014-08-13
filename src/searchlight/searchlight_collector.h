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
 * @file searchlight_collector.h
 *
 * This files describes a validated solutions collector, which is basically
 * a composition of an or-tools solution collector and SearchlighTask. Since
 * the collector can be instantiated only at the validator, this class
 * provides the necessary means to do so.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_SEARCHLIGHT_COLLECTOR_H_
#define SEARCHLIGHT_SEARCHLIGHT_COLLECTOR_H_

#include "ortools_inc.h"

namespace searchlight {

class SearchlightTask;

/**
 * This class collects tracks exacts solutions, stringifies them and them
 * to the searchlight task. It also signals to the task when the search ends.
 */
class TaskSolutionCollector : public SolutionCollector {
public:
    /**
     * Creates a new task collector.
     *
     * @param task the task to collect for
     * @param s the solver from searchlight
     */
    TaskSolutionCollector(SearchlightTask &task, Solver *s);

    /**
     * Destructor.
     */
    virtual ~TaskSolutionCollector() {}

    /**
     * A callback called at a solution.
     *
     * @return true, if we want to continue; false otherwise
     */
    virtual bool AtSolution();

    /**
     * A callback that is called when the search is finished.
     */
    virtual void ExitSearch();

    /**
     * Returns a string representation for debugging.
     *
     * @return a string representation
     */
    virtual std::string DebugString() const {
        if (prototype_.get() == NULL) {
            return "TaskSolutionCollector()";
        } else {
            return "TaskSolutionCollector(" + prototype_->DebugString() + ")";
        }
    }

    /**
     * Converts the specified solution to a string representation.
     *
     * This function basically add variable names to the values. The collector
     * is one of the few high-level classes that can do that, since it has
     * knowledge about decision variables. We assume the values correspond to
     * the variables one-to-one.
     *
     * @param vals variable values
     * @return string representation of the solution
     */
    std::string SolutionToString(const std::vector<int64_t> &vals) const;

private:
    // The searchlight engine
    SearchlightTask &task_;
};

/**
 * This is a simple wrap around the TaskSolutionCollector to make delayed
 * initialization possible. Before the collector can be used, the user
 * should call InitCollector() and then register the collector with the
 * solver after retrieving it with GetCollector().
 *
 * One use for delayed initialization is to be able to register it at the
 * validator.
 *
 */
class SearchlightCollector {
public:
    /**
     * Creates a new searchlight collector.
     *
     * @param task the main SearchlightTask
     */
    SearchlightCollector(SearchlightTask &task) :
        collector_(nullptr),
        task_(task) {}

    /**
     * Destructor.
     */
    ~SearchlightCollector() {
        delete collector_;
    }

    /**
     * Initializes this collector for the given solver.
     *
     * @param s the solver to create the collector for
     */
    void InitCollector(Solver *s) {
        collector_ = new TaskSolutionCollector(task_, s);
    }

    /**
     * Returns the collector to use with a solver.
     *
     * @return the collector to use with a solver
     */
    SolutionCollector *GetCollector() {
        return collector_;
    }

    /**
     * Returns a string representation of a solution
     * @param vals
     * @return
     */
    std::string SolutionToString(const std::vector<int64_t> &vals) const {
        assert(collector_);
        return collector_->SolutionToString(vals);
    }

private:
    // Solution collector, which will be used with the solver
    TaskSolutionCollector *collector_;

    // The main task to give the results to
    SearchlightTask &task_;
};
}
#endif /* SEARCHLIGHT_SEARCHLIGHT_COLLECTOR_H_ */

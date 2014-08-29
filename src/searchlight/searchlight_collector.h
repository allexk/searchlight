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
class SearchlightSolutionCollector : public SolutionCollector {
public:
    /**
     * Creates a new task collector.
     *
     * @param task the task to collect for
     * @param s the solver from searchlight
     */
    SearchlightSolutionCollector(SearchlightTask &task, Solver *s);

    /**
     * A callback called at a solution.
     *
     * @return true, if we want to continue; false otherwise
     */
    virtual bool AtSolution();

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
}
#endif /* SEARCHLIGHT_SEARCHLIGHT_COLLECTOR_H_ */

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
 * @file validator.h
 * This is the validator that collects possible solutions (assignments)
 * and checks them later for the exactness. The file also contains a couple
 * of utility classes used by the validator.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_VALIDATOR_H_
#define SEARCHLIGHT_VALIDATOR_H_

#include "scidb_inc.h"
#include "ortools_inc.h"

#include <condition_variable>
#include <mutex>

namespace searchlight {

class Searchlight;
class SearchlightCollector;

/**
 * This class allows users to collect solutions (assignments) and check them
 * later for exactness when some criteria are met. The validator uses a
 * separate solver, which is a duplicate of the original, to achieve this.
 */
class Validator {
public:
    /**
     * Creates a validator for the given vars of the given solver.
     *
     * Basically, it creates a duplicate solver that would periodically
     * set the vars to the validated values.
     *
     * @param solver the main search solver
     * @param var_names the names of the variables
     * @param collector the searchlight collector to grab validated solutions
     */
    Validator(const Searchlight &sl, const StringVector &var_names,
            SearchlightCollector &sl_collector);

    /**
     * Destructor.
     */
    ~Validator() {}

    /**
     * Adds a solution (assignment) to validate later.
     *
     * @param sol the solution to validate
     */
    void AddSolution(const Assignment &sol);

    /**
     * Blocks until the validor clears the validator queue. This allows the
     * caller to sync with the validator and make sure the validator can
     * accept further requests.
     */
    void Synchronize() const;

    /**
     * Signal to the validator that the search ended.
     *
     * The validator is not necessarily terminated
     * upon the exit from this function. It is assumed that this function
     * is called from another thread. In this case, the validator thread
     * should be joined to ensure proper termination.
     */
    void SignalEnd() {
        to_validate_mtx_.lock();
        search_ended_ = true;
        to_validate_mtx_.unlock();
        validate_cond_.notify_one();
    }

    /**
     * This operator calls the validator by initiating its loop. The validator
     * terminates only after somebody calls Terminate() function.
     */
    void operator()();

    /**
     * Returns the resulting status of the validator solver.
     *
     * NOTE: no thread-safe. Should be called only after the operator() call
     * completely finished, i.e., after the corresponding thread joined.
     *
     * @return the resulting status of the validator solver
     */
    bool GetValidatorSolverResult() const {
        return solver_status_;
    }

private:
    /*
     * We define the DB as a friend to grab the next portion of assignments
     */
    friend class RestoreAssignmentBuilder;

    // Returns the next portion of Assignments to validate
    AssignmentPtrVector *GetNextAssignments();

    // Check if the searclight is terminating
    bool CheckTerminate() const {
        std::unique_lock<std::mutex> validate_lock(to_validate_mtx_);
        return sl_.CheckTerminate();
    }

    // Searchlight instance
    const Searchlight &sl_;

    // Pending validations
    AssignmentPtrVector to_validate_;

    // The duplicate solver for validation
    Solver solver_;

    // The array access adapter
    AdapterPtr adapter_;

    // The prototype assignment for search variables
    Assignment search_vars_prototype_;

    // User's solution collector (don't need to delete)
    SolutionCollector *collector_;

    // Condition var to wait for solutions to validate
    mutable std::condition_variable validate_cond_;

    // Mutex to guard the validation array
    mutable std::mutex to_validate_mtx_;

    // Has the search ended
    bool search_ended_;

    // The resulting status of the validator solver
    bool solver_status_;

    /*
     * The maximum number of pending validations. If the number is exceeded the
     * main solver will block until the validator catches up.
     */
    int max_pending_validations_;

    /*
     * This parameter specifies the maximum number of assignments to check
     * before the validator makes a restart. This can be seen as periodic
     * garbage collecting, since restarting destroys elements Alloced with
     * the solver.
     */
    int restart_period_;
};


/**
 * This class is a visitor that collects all integer variables. It creates
 * a map that variable names to the addresses.
 */
class VariableFinder : public ModelVisitor {
public:
    /**
     * A var_name --> var_address map
     */
    typedef std::map<std::string, const IntVar *> StringVarMap;

    /**
     * Visits all "simple" integer variables.
     *
     * Even if the variable is a cast one, we collect its name and do not go to
     * the expression.
     *
     * @param variable the variable visited
     * @param delegate the expression this variable is casted for
     */
    virtual void VisitIntegerVariable(const IntVar* const variable,
                                      IntExpr* const delegate) override {
        if (variable->HasName()) {
            var_map_[variable->name()] = variable;
        }
    }

    /**
     * Visits "optimized" variables, like x+c, x*c, etc. and some others.
     *
     * @param variable the variable itself
     * @param operation operation name (e.g., sum, product, etc.)
     * @param value the value if any (e.g., c in x+c)
     * @param delegate the original variable (e.g., x in x+c)
     */
    virtual void VisitIntegerVariable(const IntVar* const variable,
                                      const std::string& operation, int64 value,
                                      IntVar* const delegate) override {
        /*
         * We are not interested in the delegate. We are going to work
         * with the optimized variable.
         */
        if (variable->HasName()) {
            var_map_[variable->name()] = variable;
        }
    }

    /**
     * Return a mapping of var names to var addresses.
     *
     * @return the name->address variables map
     */
    const StringVarMap &GetVarMap() const {
        return var_map_;
    }

private:
    // The map containing vars found
    StringVarMap var_map_;
};

} /* namespace searchlight */
#endif /* SEARCHLIGHT_VALIDATOR_H_ */

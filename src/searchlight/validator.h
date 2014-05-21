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

#include <boost/thread.hpp>

namespace searchlight {

class Searchlight;

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
     * @param collector the solution collector to grab validated solutions
     */
    Validator(const Searchlight &sl, const StringVector &var_names,
            SolutionCollector &collector);

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
     * Terminates the validator.
     *
     * The validator is not necesarily terminated
     * upon the exit from this function. It is assumed that this function
     * is called from another thread. In this case, the validator thread
     * should be joined to ensure proper termination.
     */
    void Terminate() {
        boost::unique_lock<boost::mutex> validate_lock(to_validate_mtx_);
        terminate_ = true;
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

    /*
     * This is a builder for UDF functions (needed to restore the model from
     * a protobuf buffer)
     */
    IntExpr* UDFBuilder(std::string name, CPModelLoader* const builder,
                          const CPIntegerExpressionProto& proto);

    // Registers builder for UDFs to restrore them from the protobuf model
    void RegisterUDFBuilder();

    // Returns the next portion of Assignments to validate
    AssignmentPtrVector *GetNextAssignments();

    // Check if the validator is terminating
    bool CheckTerminate() const {
        boost::unique_lock<boost::mutex> validate_lock(to_validate_mtx_);
        return terminate_;
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

    // User's sollution collector
    SolutionCollector &collector_;

    // Condition var to wait for solutions to validate
    boost::condition_variable validate_cond_;

    // Mutex to guard the validation array
    mutable boost::mutex to_validate_mtx_;

    // Should we terminate?
    bool terminate_;

    // The resulting status of the validator solver
    bool solver_status_;
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
     * Visits all "simple" integer variables. We collect only non-cast
     * variables that have names.
     *
     * @param variable the variable visited
     * @param delegate the expression this variable is casted for
     */
    virtual void VisitIntegerVariable(const IntVar* const variable,
                                      const IntExpr* const delegate) {
        if (!delegate && variable->HasName()) { // Ignore cast vars
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

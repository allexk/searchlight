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
     * @param sl_task current searchlight task
     * @param var_names the names of the variables
     */
    Validator(Searchlight &sl, SearchlightTask &sl_task,
            const StringVector &var_names);

    /**
     * Adds a solution (assignment) to validate later.
     *
     * @param sol the solution to validate
     */
    void AddSolution(const Assignment &sol);

    /**
     * Submits forwarded candidates for validation.
     *
     * @param cands candidates to validate
     * @param src source validator
     * @param forw_id remote forward id of the first candidate
     */
    void AddRemoteCandidates(LiteAssignmentVector &cands, InstanceID src,
            int forw_id);

    /**
     * Handles the result of the forwarder validation.
     *
     * @param id forward id
     * @param result true, if the candidate is valid; false, otherwise
     */
    void HandleForwardResult(int id, bool result);

    /**
     * Blocks until the validor clears the validator queue. This allows the
     * caller to sync with the validator and make sure the validator can
     * accept further requests.
     */
    void Synchronize() const;

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

    /**
     * Wake-ups validator if it's idle.
     *
     * Waking up validator will ensure that it can re-check the SL status and
     * terminate properly.
     */
    void WakeupIfIdle() const {
        validate_cond_.notify_one();
    }

private:
    // Info for a candidate assignment pending validation
    struct CandidateAssignment {
        // The assignment itself
        LiteVarAssignment var_asgn_;

        // Id >= 0, if the candidate is a remote; -1, local
        int forw_id_;
    };
    using CandidateVector = std::vector<CandidateAssignment>;

    /*
     * We define the DB as a friend to grab the next portion of assignments
     */
    friend class RestoreAssignmentBuilder;

    // Returns the next portion of Assignments to validate
    CandidateVector *GetNextAssignments();

    // Check if the searchlight is terminating
    bool CheckTerminate() const {
        return sl_.CheckTerminate();
    }

    /*
     *  Returns true if this validator finished locally:
     *    1) No local candidate solutions
     *    2) No outstanding forwards
     */
    bool FinishedLocally() const {
        /*
         *  No mutex here since it's called from the main validator loop,
         *  which takes care of that.
         */
        return to_validate_.empty() && forwarded_candidates_.empty();
    }

    // Sends back the result of the forward
    void SendForwardResult(int forw_id, bool result);

    // Checks if we want to forward (returns true) and forwards if we can.
    bool CheckForward(const Assignment *asgn);

    // Searchlight instance
    Searchlight &sl_;

    // Searchlight task
    SearchlightTask &sl_task_;

    // Pending validations
    CandidateVector to_validate_;

    // Info about remote candidates (local id -> (instance, remote id))
    std::unordered_map<int, std::pair<InstanceID, int>> remote_candidates_;

    // Local ids for remote candidates
    int remote_cand_id_ = 0;

    // Contains information about forwarded candidates (local_id -> solution)
    std::unordered_map<int, LiteVarAssignment> forwarded_candidates_;

    // Forward ids for candidates
    int forw_id_ = 0;

    // True, if for candidate forwarding is enabled
    bool balancing_enabled_;

    // The duplicate solver for validation
    Solver solver_;

    // The array access adapter
    AdapterPtr adapter_;

    // The prototype assignment for search variables
    Assignment search_vars_prototype_;

    // Condition var to wait for solutions to validate
    mutable std::condition_variable validate_cond_;

    // Mutex to guard the validation array
    mutable std::mutex to_validate_mtx_;

    // The resulting status of the validator solver
    bool solver_status_ = false;

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

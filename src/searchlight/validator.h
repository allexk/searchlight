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
#include "ortools_model.h"

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
     * @param coords relevant coordinates to determine forwards zones
     */
    void AddRemoteCandidates(LiteAssignmentVector &cands,
            const std::vector<std::vector<Coordinates1D>> &coords,
            InstanceID src, int forw_id);

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
        // Avoid the case where we signal before the validator starts waiting
        std::lock_guard<std::mutex> lock{to_validate_mtx_};
        validate_cond_.notify_all();
    }

    /**
     * Add a new helper for the validator.
     *
     * The request is non-binding, although the thread is supposed to
     * be reserved by the caller. The validator might decide not to add
     * any helpers, in which case the thread is returned back.
     *
     * The function will return the persistence of the added worker (if
     * it has been added). A persistent worker is supposed to be long-term, so
     * the caller might plan accordingly. A non-persistent worker checks in
     * with the task regularly to see if it needs to quit.
     *
     * @param persistent pointer to put the persistence of the worker to
     * @return true, if the helper's been started; false, otherwise
     */
    bool AddValidatorHelper(bool *persistent);

    /**
     * Update the number of candidates statistics for the specified validator.
     *
     * @param cands the number of candidates
     * @param val validator ordinal number
     */
    void UpdateCandsInfo(size_t cands, size_t val) {
        // No locking, this is just stats
        validators_cands_info_[val] = cands;
    }

private:
    // Type of validator forwarding
    enum class Forwarding {
        // No forwarding
        NONE,
        // Each active validator takes a stripe
        STRIPES,
        // Based on the chunks each validator has (static + dynamic)
        DYNAMIC
    };

    // Info for a candidate assignment pending validation
    struct CandidateAssignment {
        // The assignment itself
        LiteVarAssignment var_asgn_;

        // Id >= 0, the candidate is remote; -1, needs simulation; -2, local
        int forw_id_;
    };
    using CandidateVector = std::vector<CandidateAssignment>;

    /*
     * We define the DB as a friend to grab the next portion of assignments
     */
    friend class RestoreAssignmentBuilder;

    /**
     * ValidatorHelper is an advanced solver for validating candidates.
     *
     * The basic idea behind a helper is to off-load a bunch of candidates
     * to it to validate in parallel. It uses the main validator to
     * perform all important stuff, like forwarding and reporting.
     */
    class ValidatorHelper {
    public:
        /**
         * Unique pointer for ValidatorHelper.
         */
        using UniquePtr = std::unique_ptr<ValidatorHelper>;

        /**
         * Constructs a new validator helper.
         *
         * @param id unique id if the helper
         * @param parent main validator
         * @param init_assignment assignment to init the vars
         */
        ValidatorHelper(int id, Validator &parent,
                const Assignment &init_assignment);

        /**
         * Destructor.
         *
         * If the thread is still joinable (validator joins them lazily) then
         * it will be joined.
         */
        ~ValidatorHelper() {
            Prepare();
        }

        /**
         * Main loop for the helper's solver.
         */
        void operator()();

        /**
         * Starts the specified workload in a separate thread.
         *
         * The helper must not have the current workload and it must not
         * be running in a thread.
         *
         * @param workload workload to validate
         */
        void RunWorkload(CandidateVector &&workload);

        /**
         * Set the persistence of the helper.
         *
         * A persistent helper lasts until it reaches below the low watermark.
         * Then it checks in with the task. A non-persistent helper checks in
         * every time when below the high watermark.
         *
         * All helpers are non-persistent by default.
         *
         * @param pers persistance value
         */
        void SetPersistent(bool pers) {
            persistent_ = pers;
        }

    private:
        /*
         * Internal helper for setting up the workload.
         */
        void SetWorkload(CandidateVector &&workload) {
            workload_ = std::move(workload);
        }

        /*
         * Internal helper to prepare the helper for use.
         */
        void Prepare() {
            if (thr_.joinable()) {
                thr_.join();
            }
        }

        // Parent validator
        Validator &parent_;

        // Validation workload
        CandidateVector workload_;

        // Solver
        Solver solver_;

        // Adapter to use for UDFs
        AdapterPtr adapter_;

        // The prototype assignment for search variables
        Assignment prototype_;

        // Thread to run the helper
        std::thread thr_;

        // Persistent helper (lasts until the low watermark)
        bool persistent_ = false;

        // Id of the helper
        const int id_;
    };

    // Returns the next portion of Assignments to validate
    CandidateVector *GetNextAssignments();

    // Check if the searchlight is terminating
    bool CheckTerminate() const {
        return sl_.CheckTerminate();
    }

    // Tries to dispatch a helper; nullptr -- not possible
    ValidatorHelper *DispatchValidatorHelper();

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
        return to_validate_total_ == 0 && forwarded_candidates_.empty() &&
                free_validator_helpers_.size() == validator_helpers_.size();
    }

    // Sends back the result of the forward
    void SendForwardResult(int forw_id, bool result);

    // Checks if we want to forward (returns true) and forwards if we can.
    bool CheckForward(const CoordinateSet &chunks, const Assignment *asgn);

    // Determines the assignment's zone and pushes the candidate there
    void PushToLocalZone(const CoordinateSet &chunks, const Assignment *asgn);

    // Determine the local zone to use
    template <typename CoordinatesSequence>
    size_t DetermineLocalZone(const CoordinatesSequence &chunks) const;

    // Pushes candidate to to_validate_
    void PushCandidate(CandidateAssignment &&asgn, size_t zone);

    // Returns a new assignment workload
    CandidateVector GetWorkload();

    // Adds a helper to the pool (assumes the mutex is locked)
    bool AddValidatorHelperInt(bool *persistent);

    // Return the number of local zones
    size_t LocalZonesNumber() const {
        // one "zone" for new, non-simulated, candidates
        return to_validate_.size() - 1;
    }

    // Searchlight instance
    Searchlight &sl_;

    // Searchlight task
    SearchlightTask &sl_task_;

    // Validator model for later cloning
    CPModelProto validator_model_;

    // Pending validations and their count (vector of zones)
    std::vector<std::deque<CandidateVector>> to_validate_;
    size_t to_validate_total_ = 0;
    // Zones in the MRU order (we use vector since that might be faster even
    // for in-the-middle erases than list splicing.
    // The MRU zone goes last!
    std::vector<int> zones_mru_;

    // Periodicity of sending the number of candidates updates and the info
    std::size_t send_info_period_;
    std::size_t last_sent_cands_num_ = 0;
    std::vector<size_t> validators_cands_info_;

    // Info about remote candidates (local id -> (instance, remote id))
    std::unordered_map<int, std::pair<InstanceID, int>> remote_candidates_;

    // Local ids for remote candidates
    int remote_cand_id_ = 0;

    // Contains information about forwarded candidates (local_id -> solution)
    std::unordered_map<int, LiteVarAssignment> forwarded_candidates_;

    // Forward ids for candidates
    int forw_id_ = 0;

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

    /*
     * Type of forwarding used.
     */
    Forwarding forw_type_;

    // Logical id of this validator among active ones (cached for perfomance)
    int my_logical_id_;

    // Initialized validator helpers
    std::unordered_map<int, ValidatorHelper::UniquePtr> validator_helpers_;

    // Free validator helpers (finished their workloads)
    std::deque<int> free_validator_helpers_;

    // Maximum number of helpers allowed
    size_t max_helpers_allowed_;

    // Number of assignments to off-load to a helper
    size_t helper_workload_;

    // Low- and high-watermark to determine the necessity of help
    size_t low_watermark_, high_watermark_;

    // Do we use dynamic helper scheduling?
    bool dynamic_helper_scheduling_;

    // Chunk zones to determine forwards and local candidate batches
    SearchArrayDesc::ChunkZones chunk_zones_;
};


/**
 * This class is a visitor that collects all integer variables. It creates
 * a map that variable names to the addresses.
 */
class VariableFinder : public ModelVisitor {
public:
    /**
     * A var_name --> var_address map.
     *
     * The IntVar is const since that is the way they are traversed at model
     * visitors.
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

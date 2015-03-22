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
 * @file searchlight_task.h
 *
 * This file contains class definitions governing the search. It includes
 * the array for delivering the result upstream and the searchlight task,
 * which runs the search in am individual thread. Also, it contains a
 * monitor that saves confirmed results into an array's buffer.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_TASK_H_
#define SEARCHLIGHT_TASK_H_

#include "ortools_inc.h"
#include "scidb_inc.h"
#include "searchlight.h"
#include "searchlight_messenger.h"

#include <mutex>
#include <condition_variable>
#include <functional>
#include <unordered_set>

namespace searchlight {

/**
 * This class incapsulates the searchlight process and a monitor that collects
 * results and put them into a queue. This class is also repsonsible for
 * resolving and establishing the task. It is supposed to be running in a
 * thread. Notice, that the searchlight's validator will run another thread,
 * so its going to be two (2) threads per task.
 */
class SearchlightTask {
public:
    /**
     * Constructs a searchlight task.
     *
     * @param library_name the name of the DLL library with the task
     *  (without .so extension)
     * @param task_name the name of the task (basically, the function name)
     * @param config_file_name config file for the task (parsed by searchlight)
     * @param data data array
     * @param samples vector of sample arrays
     * @param query current query
     */
    SearchlightTask(const std::string &library_name,
            const std::string &task_name, const std::string &config_file_name,
            ArrayPtr &data, const ArrayPtrVector &samples,
            const boost::shared_ptr<Query> &query);

    /**
     * Returns the next solution if any. Note, it might block for a
     * considrable period of time waiting for a solution.
     *
     * @return the next solution or an empty string if there are no left
     */
    std::string GetNextSolution();

    /**
     * Terminate the current task. In a common (successful) case the solver
     * and validator will already be committed, so no action would be
     * required. In a bad case (i.e., an error or cancelled query) we might
     * need to stop both the solver and validator. In this case we just set
     * the status, and they will pick it up from there.
     */
    void Terminate() {
        if (InstanceActive(my_instance_id_) &&
                searchlight_.GetStatus() != Searchlight::Status::COMMITTED) {
            searchlight_.Terminate();
        }
    }

    /**
     * Handle error in one of the Searchlight components.
     *
     * Basically, that causes propagation of the error to other
     * components via the standard Terminate mechanism and termination
     * of the task itself.
     *
     * @param error error causing termination
     */
    void HandleSearchlightError(const scidb::Exception::Pointer &error) {
        std::lock_guard<std::mutex> lock{mtx_};
        sl_error_ = error;
        Terminate();
        queue_cond_.notify_one();
    }

    /**
     * Determines if the task terminated erroneously.
     * @return true, if terminated on error; false, otherwise
     */
    bool ErrorTerminate() const {
        return sl_error_;
    }

    /**
     * Returns this task's instance id.
     *
     * @return this task's instance id
     */
    InstanceID GetInstanceID() const {
        return my_instance_id_;
    }

    /**
     * Returns the total number of instances participating in the query.
     *
     * @return total number of query instances
     */
    int GetQueryInstanceCount() const {
        return query_instance_count_;
    }

    /**
     * Reports an idle solver to the coordinator.
     *
     * If the solver is local to the coordinator, the coordinator will
     * handle it locally, no messages involved.
     *
     * If the <code>postponed</code> parameter is true, then the solver's
     * reporting might be postponed due to dynamic scheduling, if the
     * validator need more helpers. If the dynamic scheduling is off, the
     * parameter doesn't make any effect.
     *
     * @param id solver id
     * @param deferable true, if the solver's reporting can be postponed
     */
    void ReportIdleSolver(uint64_t solver_id, bool deferable);

    /**
     * Reports locally finished validator.
     *
     * Validator is finished locally if it does not have any local requests
     * coming and no outstanding forwards waiting on.
     */
    void ReportFinValidator();

    /**
     * Rejects help from other solvers, returning helpers to the coordinator.
     *
     * @param helpers helpers to reject
     * @param solver_id id of the solver rejecting help
     * @param hard true, if hard reject; false, if soft
     */
    void RejectHelp(const std::vector<uint64_t> &helpers, uint64_t solver_id,
            bool hard);

    /**
     * Dispatches work to another solver.
     *
     * This function is used by a solver to off-load some of its work to
     * another solver that was dispatched by the coordinator as a helper.
     *
     * As one of its steps it notifies the coordinator that the helper
     * has been accepted.
     *
     * @param work assignments to send
     * @param dest_solver id of the solver to send the work to
     */
    void DispatchWork(const LiteAssignmentVector &work,
            uint64_t dest_solver);

    /**
     * Forwards candidates to another validator.
     *
     * The destination validator will reply with the same id when it's done.
     *
     * @param cands candidates to forward
     * @param zones candidate zones
     * @param dest destination validator
     * @param forw_id id of the forward
     */
    void ForwardCandidates(const LiteAssignmentVector &cands,
            const std::vector<int64_t> &zones,
            InstanceID dest, uint64_t forw_id) const;

    /**
     * Sends result of the balancing.
     *
     * @param dest destination Searchlight
     * @param id id of the balancing load
     * @param result result status of the load
     */
    void SendBalanceResult(InstanceID dest, uint64_t id, bool result) const;

    /**
     * Broadcast validator's info.
     *
     * Currently we broadcast only the number of candidates in the queue.
     *
     * @param cands_num the number of candidates
     */
    void BroadcastValidatorInfo(size_t cands_num) const;

    /**
     * Returns query context of this task.
     *
     * @return task's query context
     */
    boost::shared_ptr<Query> GetQueryContext() const {
        return Query::getValidQueryPtr(query_);
    }

    /**
     * Reports a real (not candidate) solution.
     *
     * This function will either cause the solution to be reported to the
     * user (via SciDb) or be sent to the coordinator.
     *
     * @param values solution values
     */
    void ReportSolution(const std::vector<int64_t> &values);

    /**
     * Returns the property tree containing configuration options.
     *
     * @return property tree with the config
     */
    const SearchlightConfig &GetConfig() const {
        return config_;
    }

    /**
     * Returns a list of instances containing active validators.
     *
     * Active validators are specified in the configuration file.
     *
     * @return list of active validator instances
     */
    const std::vector<InstanceID> &GetActiveValidators() const {
        return active_validators_;
    }

    /**
     * Return number of solvers at the instance.
     *
     * @param id instance id
     * @return number of solvers at the instance
     */
    int GetSolverNum(InstanceID id) const {
        auto it = std::find(active_solver_instances_.begin(),
                        active_solver_instances_.end(), id);
        if (it != active_solver_instances_.end()) {
            return active_solver_num_[it - active_solver_instances_.begin()];
        } else {
            return 0;
        }
    }

    /**
     * Checks if the solver is active at the instance.
     *
     * @param id instance to check
     * @return true, if the solver is active; false, otherwise
     */
    bool SolverActive(InstanceID id) const {
        return std::find(active_solver_instances_.begin(),
                active_solver_instances_.end(), id) !=
                        active_solver_instances_.end();
    }

    /**
     * Checks if the specified instance is active.
     *
     * An instance is active if it has an active solver or validator.
     *
     * @param id instance to check
     * @return true, if the instance is active; false, otherwise
     */
    bool InstanceActive(InstanceID id) const {
        const bool solver_active = SolverActive(id);
        const bool validator_active = std::find(active_validators_.begin(),
                active_validators_.end(), id) != active_validators_.end();
        return solver_active || validator_active;
    }

    /**
     * Returns total solvers number.
     * @return
     */
    size_t GetTotalSolversNum() const {
        size_t total = 0;
        for (auto s: active_solver_num_) {
            total += s;
        }
        return total;
    }

    /**
     * Return global orinal solver number.
     *
     * @param id solver's id
     * @return global ordinal id for the solver
     */
    size_t GetGlobalOrdinalSolverId(uint64_t id) const;

    /**
     * Start the search process.
     *
     * We start search only if we have active solvers at this instance. If
     * there is an active validator, it has already been started by the
     * prepare function in Searchlight.
     */
    void StartSearch() {
        if (SolverActive(my_instance_id_)) {
            searchlight_.StartSolvers();
        }
    }

    /**
     * Reserve a thread if one is available.
     *
     * It is assumed that the caller respects the return value. If it is
     * false, that means no available threads are available. True means
     * a new thread can be created and given work to do.
     *
     * @return true, if the reservation was a success; false, otherwise
     */
    bool ReserveThread() {
        if (!dynamic_scheduling_) {
            return true;
        }

        std::lock_guard<std::mutex> lock{mtx_};
        if (threads_available_ > 0) {
            threads_available_--;
            return true;
        } else {
            return false;
        }
    }

    /**
     * Free one of the reserved threads.
     */
    void FreeThread();

    /**
     * Check if there are somesolver jobs pending and some solvers are
     * waiting to be reported as idle.
     *
     * If the function detects pending jobs, that means the caller should
     * free their own thread ASAP. While this is non-binding, it should be
     * respected.
     *
     * The additional check_idle_solvers parameter signals the task to
     * check if some of the solvers are waiting to be reported as idle. This
     * is primarily useful for persistent helpers, who want to signal they
     * no longer needed for flood control. If there are indeed idle solvers,
     * one of them will be reported to the coordinator.
     *
     * @return true, if there are solver jobs; false, otherwise
     */
    bool PendingSolverJobs(bool check_idle_solvers);

private:
    // Contains information about the state of distributed search.
    struct DistributedSearchInfo {
        // Currently busy solvers
        std::unordered_set<uint64_t> busy_solvers_;

        // Currently idle solvers
        std::unordered_set<uint64_t> idle_solvers_;

        // Solvers needing help
        struct {
            // LRU queue of solvers that need help
            std::list<uint64_t> help_reqs_;

            // Index to quickly find and remove instances from the list
            std::unordered_map<uint64_t, std::list<uint64_t>::iterator> index_;

            // Erases solver id from helpees
            void Erase(uint64_t id) {
                const auto iter = index_.find(id);
                if (iter != index_.end()) {
                    help_reqs_.erase(iter->second);
                    index_.erase(iter);
                }
            }

            // Add solver id to helpees
            void Add(uint64_t id) {
                // Might be already here if accept ack came too late.
                // In this case we might've been added by another helpee...
                if (index_.find(id) == index_.end()) {
                    help_reqs_.push_back(id);
                    index_.emplace(id, --help_reqs_.end());
                }
            }
        } helpees_;

        // Number of validators still busy with the job
        int busy_validators_count_;
    };

    // Info about distributed search (non-NULL only at the coordinator)
    std::unique_ptr<DistributedSearchInfo> distr_search_info_;

    // Returns the next solver that requires help
    uint64_t GetHelpee();

    // Resolves the task in the DLL
    void ResolveTask(const std::string &lib_name, const std::string &task_name);

    // Handles an idle solver and possible dispatches it to another solver
    void HandleIdleSolver(uint64_t id);

    // Handles a locally finished validator
    void HandleFinValidator(InstanceID id);

    // Handles the end-of-search messsage
    void HandleEndOfSearch();

    // Handles the commit message
    void HandleCommit();

    // Are we expecting more results?
    bool ExpectingResults() const {
        assert(distr_search_info_);
        return distr_search_info_->busy_validators_count_ > 0;
    }

    // Adds solution from a remote instance into the queue
    void HandleRemoteSolution(InstanceID inst,
            const google::protobuf::Message *msg);

    // Handles control messages
    void HandleControlMessage(InstanceID inst,
            const google::protobuf::Message *msg);

    // Handles balance messages
    void HandleBalanceMessage(InstanceID inst,
            const google::protobuf::Message *msg);

    // Handles remote workload for the solver (helper)
    void HandleRemoteLoad(const SearchlightBalance &msg, uint64_t helper);

    // Handles forwarded candidates
    void HandleForwards(const SearchlightBalance &msg, InstanceID src);

    // Handles rejected help
    void HandleRejectHelp(uint64_t src,
            const std::vector<uint64_t> &helpers, bool hard);

    // Broadcasts control message to finish the main search
    void BroadcastFinishSearch();

    // Broadcasts control message to commit Searchlight
    void BroadcastCommit();

    // Checks if somebody needs help and dispatches a helper
    void CheckForHelpees();

    // Dispatches helper to another solver
    void DispatchHelper(uint64_t helper, uint64_t dest);

    // Handles accept help message
    void HandleAcceptHelper(uint64_t helper);

    // Handles a new helper (note: it will be rejected/accepted later)
    void HandleHelper(uint64_t helper, uint64_t helpee);

    // Reads config from the specified file. Only JSON is supported for now.
    void ReadConfig(const std::string &file_name);

    // Handler when a new thread is available
    void HandleFreeThread();

    // Dispatches a new load to a solver or stores it until a thread is ready.
    // Note: the load will be unusable after the call
    void DispatchOrStoreLoad(LiteAssignmentVector &load,
            uint64_t dest_solver);

    // The main DLL Handler (we need it to be destroyed last!)
    DLLHandler dll_handler_;

    // Number of instances participating in the query
    int query_instance_count_;

    // Our instance id
    InstanceID my_instance_id_;

    // Threads available for running
    size_t threads_available_ = 0;

    // Are we using dynamic scheduling for validator helpers/solvers?
    bool dynamic_scheduling_;

    // The main sl instance
    Searchlight searchlight_;

    // Task config
    SearchlightConfig config_;

    /*
     * Active solver instances, solvers per instance and validators.
     *
     * Note, there might be multiple sovers per instance, but only a
     * single validator.
     */
    std::vector<InstanceID> active_solver_instances_;
    std::vector<int> active_solver_num_;
    std::vector<InstanceID> active_validators_;

    // The task function
    SLTaskFunc task_;

    // Stringified solutions.
    std::list<std::string> solutions_queue_;

    // Set of solvers waiting to be reported idle (due to dynamic scheduling)
    std::unordered_set<uint64_t> pending_idle_solvers_;

    // Pending external loads waiting for a solver
    std::unordered_map<uint64_t, LiteAssignmentVector> pending_assgns_;

    // Mutex to protect shared task structures
    std::mutex mtx_;

    // Condition to facilitate waiting for solutions
    std::condition_variable queue_cond_;

    // Error in the searchlight thread, if any
    scidb::Exception::Pointer sl_error_;

    // Current query
    const boost::weak_ptr<Query> query_;
};

/**
 *  A shared pointer for the SearclightTask.
 */
typedef std::shared_ptr<SearchlightTask> SearchlightTaskPtr;

/**
 * This class represents an array for returning SL results upstream.
 * This is a stream array, which means it is single-pass. The latter is
 * because we do not know the number of results beforehand and cannot
 * let the user to randomly walk through it. For the same reason its scheme's
 * only dimension is unbounded on the right.
 *
 * This array also handles the searchlight task by running it in a
 * separate thread and terminating when necessary. The search starts running
 * when the first result is requested.
 *
 * Since the array is stream-based, the results can be delivered in an online
 * fashion, which (hopefully) is what SciDb should be doing, since
 * materialization is unnecessary.
 */
class SearchlightResultsArray : public StreamArray {
public:
    /**
     * Creates a new resulting array.
     * @param sl_task the searchlight task
     * @param desc the descriptor of the array
     * @param query the query context of this array
     */
    SearchlightResultsArray(SearchlightTaskPtr sl_task, ArrayDesc &desc,
            const boost::shared_ptr<Query> &query) :
        StreamArray(desc, false),
        sl_task_(sl_task),
        desc_(desc),
        res_count_(0) {

        // Start search
        sl_task_->StartSearch();

        // we need the context to create the array
        assert(query);
        _query = query;
    }

    /**
     * Destructor.
     */
    virtual ~SearchlightResultsArray() {
        sl_task_->Terminate();
    }

    /**
     * Creates a descriptor for this type of arrays. Basically, it is a
     * single dimension [1, *] array with a single string attribute. The
     * dimension index equals the result's number, and the attribute
     * contains the stringified result.
     *
     * @return a suitable descruptor for the SL result array
     */
    static ArrayDesc GetArrayDesc();

protected:
    // Creates a chunk with a new result (one result per chunk)
    virtual ConstChunk const* nextChunk(AttributeID attId, MemChunk& chunk);

private:
    // The main SL task
    SearchlightTaskPtr sl_task_;

    // This array descriptor
    ArrayDesc desc_;

    // Current results count
    uint64_t res_count_;
};

} /* namespace searchlight */
#endif /* SEARCHLIGHT_TASK_H_ */

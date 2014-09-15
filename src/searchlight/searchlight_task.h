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
            const boost::shared_ptr<Query> &query) :
                query_instance_count_(query->getInstancesCount()),
                my_instance_id_(query->getInstanceID()),
                searchlight_(*this, task_name, dll_handler_),
                query_(query) {

        // Fill in distributed search info
        if (query->getCoordinatorID() == scidb::COORDINATOR_INSTANCE) {
            distr_search_info_.reset(new DistributedSearchInfo);
            for (int i = 0; i < query_instance_count_; i++) {
                distr_search_info_->helpees_.Add(i);
                distr_search_info_->busy_solvers_.insert(i);
            }
            distr_search_info_->busy_validators_count_ = query_instance_count_;
        }

        ResolveTask(library_name, task_name);
        searchlight_.ReadConfig(config_file_name);
        searchlight_.RegisterArray(data, samples);

        using std::placeholders::_1; // we have a clash with Boost
        using std::placeholders::_2; // we have a clash with Boost
        SearchlightMessenger::getInstance()->RegisterUserMessageHandler(
             query,
             SearchlightMessenger::mtSLSolution,
             std::bind(&SearchlightTask::HandleRemoteSolution, this, _1, _2));
        SearchlightMessenger::getInstance()->RegisterUserMessageHandler(
                query,
                SearchlightMessenger::mtSLControl,
                std::bind(&SearchlightTask::HandleControlMessage, this, _1, _2));
        SearchlightMessenger::getInstance()->RegisterUserMessageHandler(
                query,
                SearchlightMessenger::mtSLBalance,
                std::bind(&SearchlightTask::HandleBalanceMessage, this, _1, _2));

        // Tasks just prepare Searchlight, solver is still idle, no threads
        task_(&searchlight_);
    }

    /**
     * Operator for running the task. Basically, just calls the function from
     * the DLL.
     */
    void operator()();

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
        if (searchlight_.GetStatus() != Searchlight::Status::COMMITTED) {
            searchlight_.Terminate();
            sl_cond_.notify_one();
        }
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
     */
    void ReportIdleSolver();

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
     * @param hard true, if hard reject; false, if soft
     */
    void RejectHelp(const std::vector<InstanceID> &helpers, bool hard);

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
     * @param solver id of the solver to send the work to
     */
    void DispatchWork(const LiteAssignmentVector &work,
            InstanceID solver);

    /**
     * Forwards candidates to another validator.
     *
     * The destination validator will reply with the same id when it's done.
     *
     * @param cands candidates to forward
     * @param dest destination validator
     * @param forw_id id of the forward
     */
    void ForwardCandidates(const LiteAssignmentVector &cands,
            InstanceID dest, int forw_id) const;

    /**
     * Sends result of the balancing.
     *
     * @param dest destination Searchlight
     * @param id id of the balancing load
     * @param result result status of the load
     */
    void SendBalanceResult(InstanceID dest, int id, bool result) const;

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

private:
    // Make it a friend to modify the queue
    friend class SearchlightSolutionCollector;

    // Contains information about the state of distributed search.
    struct DistributedSearchInfo {
        // Currently busy solvers
        std::unordered_set<InstanceID> busy_solvers_;

        // Currently idle solvers
        std::unordered_set<InstanceID> idle_solvers_;

        // Solvers needing help
        struct {
            // LRU queue of solvers that need help
            std::list<InstanceID> help_reqs_;

            // Index to quickly find and remove instances from the list
            std::unordered_map<InstanceID, std::list<InstanceID>::iterator> index_;

            // Erases instance id from helpees
            void Erase(InstanceID id) {
                const auto iter = index_.find(id);
                if (iter != index_.end()) {
                    help_reqs_.erase(iter->second);
                    index_.erase(iter);
                }
            }

            // Add instance id to helpees
            void Add(InstanceID id) {
                assert(!index_.count(id));
                help_reqs_.push_back(id);
                index_.emplace(id, --help_reqs_.end());
            }
        } helpees_;

        // Number of validators still busy with the job
        int busy_validators_count_;
    };

    // Info about distributed search (non-NULL only at the coordinator)
    std::unique_ptr<DistributedSearchInfo> distr_search_info_;

    // Returns the next instance that requires help
    InstanceID GetHelpee();

    // Resolves the task in the DLL
    void ResolveTask(const std::string &lib_name, const std::string &task_name);

    // Handles an idle solver and possible dispatches it to another solver
    void HandleIdleSolver(InstanceID id);

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
    void HandleRemoteLoad(const SearchlightBalance &msg);

    // Handles forwarded candidates
    void HandleForwards(const SearchlightBalance &msg, InstanceID src);

    // Handles rejected help
    void HandleRejectHelp(InstanceID src,
            const std::vector<InstanceID> &helpers, bool hard);

    // Broadcasts control message to finish the main search
    void BroadcastFinishSearch();

    // Broadcasts control message to commit Searchlight
    void BroadcastCommit();

    // Checks if somebody needs help and dispatches a helper
    void CheckForHelpees();

    // Dispatches helper to another node
    void DispatchHelper(InstanceID helper, InstanceID dest);

    // Handles accept help message
    void HandleAcceptHelper(InstanceID helper);

    // Handles a new helper (note: it will be rejected/accepted later)
    void HandleHelper(InstanceID helper);

    // Type of the task function (called from the library)
    typedef void (*SLTaskFunc)(Searchlight *);

    // The main DLL Handler (we need it to be destroyed last!)
    DLLHandler dll_handler_;

    // Number of instances participating in the query
    int query_instance_count_;

    // Our instance id
    InstanceID my_instance_id_;

    // The main sl instance
    Searchlight searchlight_;

    // The task function
    SLTaskFunc task_;

    // Stringified solutions.
    std::list<std::string> solutions_queue_;

    // Mutex to protect shared task structures
    std::mutex mtx_;

    // Condition to facilitate waiting for solutions
    std::condition_variable queue_cond_;

    // Condition to for the solver to wait for work
    std::condition_variable sl_cond_;

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
        res_count_(0),
        sl_thread_(std::ref(*sl_task_)) {

        // we need the context to create the array
        assert(query);
        _query = query;
    }

    /**
     * Destructor.
     */
    virtual ~SearchlightResultsArray() {
        sl_task_->Terminate();
        sl_thread_.join();
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

    // The main sl thread
    std::thread sl_thread_;
};

} /* namespace searchlight */
#endif /* SEARCHLIGHT_TASK_H_ */

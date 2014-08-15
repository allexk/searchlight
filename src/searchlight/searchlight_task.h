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
#include "searchlight_collector.h"
#include "searchlight_messenger.h"

#include <mutex>
#include <condition_variable>
#include <functional>

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
                active_instance_count_(query->getInstancesCount()),
                searchlight_(task_name, active_instance_count_,
                        query->getInstanceID(), dll_handler_),
                collector_(*this),
                query_(query) {

        ResolveTask(library_name, task_name);
        searchlight_.RegisterCollector(&collector_);
        searchlight_.ReadConfig(config_file_name);
        searchlight_.RegisterArray(data, samples);

        using std::placeholders::_1; // we have a clash with Boost
        using std::placeholders::_2; // we have a clash with Boost
        SearchlightMessenger::getInstance()->RegisterUserMessageHandler(
                query,
                SearchlightMessenger::mtSLSolution,
                std::bind(&SearchlightTask::AddRemoteSolution, this, _1, _2));

        // Tasks just prepare Searchlight, solver is still idle, no threads
        task_(&searchlight_);
    }

    /**
     * Operator for running the task. Basically, just calls the function from
     * the DLL.
     */
    void operator()() {
        try {
            searchlight_.Solve();
        } catch (const scidb::Exception &ex) {
            /*
             * std::thread would call std::terminate(), so
             *  we have to take care of proper error reporting
             */
            sl_error_ = ex.copy();
            OnFinishSearch();
        } catch (const std::exception &e) {
            // Catch other C++ and library exceptions and translate them
            sl_error_ = (SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                    SCIDB_LE_ILLEGAL_OPERATION) << e.what()).copy();
            OnFinishSearch();
        }
    }

    /**
     * Returns the next solution if any. Note, it might block for a
     * considrable period of time waiting for a solution.
     *
     * @return the next solution or an empty string if there are no left
     */
    std::string GetNextSolution();

    /**
     * Terminate the current task. It passes along the requirement to the
     * SL search, which will terminate in a (short) period of time. Then, the
     * thread can be joined, unless something is holding it in the DLL's
     * task function, in which case there is no way to deal with that.
     */
    void Terminate() {
        searchlight_.Terminate();
    }

    /**
     * Determines if the task terminated erroneously.
     * @return true, if terminated on error; false, otherwise
     */
    bool ErrorTerminate() const {
        return sl_error_;
    }

private:
    // Make it a friend to modify the queue
    friend class TaskSolutionCollector;

    // Resolves the task in the DLL
    void ResolveTask(const std::string &lib_name, const std::string &task_name);

    // Either add solution to the queue or send it out
    void ReportSolution(const std::vector<int64_t> &values);

    // Signals that the search just finished
    void OnFinishSearch();

    // Are we expecting more results?
    bool ExpectingResults() const {
        return active_instance_count_ != 0;
    }

    void AddRemoteSolution(InstanceID inst,
            const google::protobuf::Message *msg);

    // Type of the task function (called from the library)
    typedef void (*SLTaskFunc)(Searchlight *);

    // The main DLL Handler (we need it to be destroyed last!)
    DLLHandler dll_handler_;

    // The number of active instances (makes sense only at the coordinator)
    int active_instance_count_;

    // The main sl instance
    Searchlight searchlight_;

    // The solution collector for exact results
    SearchlightCollector collector_;

    // The task function
    SLTaskFunc task_;

    // Stringified solutions.
    std::list<std::string> solutions_queue_;

    // Mutex to protect the solution queue
    std::mutex queue_mtx_;

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

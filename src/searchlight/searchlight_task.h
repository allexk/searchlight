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

namespace searchlight {

class SearchlightTask;

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
     * @param task_params params for the task (parsed by the task itself)
     */
    SearchlightTask(const std::string &library_name,
            const std::string &task_name, const std::string &task_params) :
                searchlight_(task_name),
                collector_(*this),
                task_params_(task_params),
                search_ended_(false) {

        ResolveTask(library_name, task_name);
        searchlight_.RegisterCollector(&collector_);
    }

    /**
     * Destructor.
     */
    ~SearchlightTask() {
        dlclose(dl_lib_handle_);
    }


    /**
     * Operator for running the task. Basically, just calls the function from
     * the DLL.
     */
    void operator()() {
        try {
            task_(&searchlight_, task_params_);
        } catch (const scidb::Exception &ex) {
            /*
             *  boost::thread would call std::terminate(), so
             *  we have to take care of proper error reporting
             */
            sl_error_ = ex.copy();
            OnFinishSearch();
        }
    }

    /**
     * Returns the searchlight engine instance.
     *
     * @return the searchlight engine instance
     */
    Searchlight &GetSearchlight() {
        return searchlight_;
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

    // Adds a new solution to the queue
    void AddSolution(const std::string &sol) {
        queue_mtx_.lock();
        solutions_queue_.push_back(sol);
        queue_mtx_.unlock();
        queue_cond_.notify_one();
    }

    // Signals that the search just finished
    void OnFinishSearch() {
        queue_mtx_.lock();
        search_ended_ = true;
        queue_mtx_.unlock();
        queue_cond_.notify_one();
    }

    // Type of the task function (called from the library)
    typedef void (*SLTaskFunc)(Searchlight *, const std::string &);

    // The main sl instance
    Searchlight searchlight_;

    // The solution collector for exact results
    SearchlightCollector collector_;

    // Task and its params
    SLTaskFunc task_;
    const std::string task_params_;

    // dll-based stuff
    void *dl_lib_handle_;

    // Stringified solutions.
    std::list<std::string> solutions_queue_;

    // Have the search ended?
    bool search_ended_;

    // Mutex to protect the solution queue
    boost::mutex queue_mtx_;

    // Condition to facilitate waiting for solutions
    boost::condition_variable queue_cond_;

    // Error in the searchlight thread, if any
    scidb::Exception::Pointer sl_error_;
};

/**
 *  A shared pointer for the SearclightTask.
 */
typedef boost::shared_ptr<SearchlightTask> SearchlightTaskPtr;

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
     */
    SearchlightResultsArray(SearchlightTaskPtr sl_task, ArrayDesc &desc) :
        StreamArray(desc, false),
        sl_task_(sl_task),
        desc_(desc),
        res_count_(0),
        sl_thread_(NULL) {}

    /**
     * Destructor.
     */
    virtual ~SearchlightResultsArray() {
        if (sl_thread_) {
            sl_task_->Terminate();
            sl_thread_->join();
            delete sl_thread_;
        }
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
    boost::thread *sl_thread_;
};

} /* namespace searchlight */
#endif /* SEARCHLIGHT_TASK_H_ */

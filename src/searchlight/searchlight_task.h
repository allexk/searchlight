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

#include "scidb_inc.h"
#include "ortools_inc.h"

namespace searchlight {

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
     */
    TaskSolutionCollector(SearchlightTask &task) :
        SolutionCollector(&(task.GetSearchlight().GetSolver())),
        task_(task) {}

    /**
     * A callback called at a solution.
     *
     * @return true, if we want to continue; false otherwise
     */
    virtual bool AtSolution();

    /**
     * A callback that is called when the search is finished.
     */
    virtual void ExitSearch() {
        task_.OnFinishSearch();
    }

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

private:
    // The searchlight engine
    SearchlightTask &task_;
};


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
                collector_(searchlight_),
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
        task_(&searchlight_, task_params_);
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
    TaskSolutionCollector collector_;

    // Task and its params
    SLTaskFunc task_;
    const std::string task_params_;

    // so-based stuff
    void *dl_lib_handle_;

    // Stringified solutions.
    std::list<std::string> solutions_queue_;

    // Have the search ended?
    bool search_ended_;

    // Mutex to protect the solution queue
    boost::mutex queue_mtx_;

    // Condition to facilitate waiting for solutions
    boost::condition_variable queue_cond_;
};

} /* namespace searchlight */
#endif /* SEARCHLIGHT_TASK_H_ */

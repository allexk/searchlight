/*
 * Copyright 2014, Brown University, Providence, RI.
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
 * @file searchlight_task.cpp
 * The implementation of the main searchlight task classes.
 *
 * @author Alexander Kalinin
 */

#include "searchlight_task.h"

#include <dlfcn.h>

namespace searchlight {

bool TaskSolutionCollector::AtSolution() {
    // we can reuse the same prototype and do not have to store the values
    const Assignment::IntContainer &vars = prototype_.get()->IntVarContainer();

    // stringify the solution
    std::ostringstream sol_string;
    for (size_t i = 0; i < vars.Size(); i++) {
        const IntVar *v = vars.Element(i).Var();
        sol_string << v->name() << "=" << v->Value();
        if (i != vars.Size() - 1) {
            sol_string << ", ";
        }
    }

    task_.AddSolution(sol_string.str());
    return true;
}

void SearchlightTask::ResolveTask(const std::string &lib_name,
        const std::string &task_name) {
    // loading the task library
    const std::string &plugins_dir = scidb::Config::getInstance()->
            getOption<std::string>(scidb::CONFIG_PLUGINS);
    std::string lib_path = plugins_dir + "/" + lib_name + ".so";
    dl_lib_handle_= dlopen(lib_path.c_str(), RTLD_LAZY | RTLD_LOCAL);
    if (!dl_lib_handle_) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Cannot load the task library: " << lib_name;
    }

    // look up the task
    task_ = (SLTaskFunc)dlsym(dl_lib_handle_, task_name.c_str());
    // We should check an error via dlerror, but NULL checking is fine
    if (!task_) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Cannot find an SL task function, name=" << task_name;
    }
}

std::string SearchlightTask::GetNextSolution() {
    boost::unique_lock<boost::mutex> lock(queue_mtx_);
    while (solutions_queue_.empty() && !search_ended_) {
        queue_cond_.wait(lock);
    }

    if (solutions_queue_.empty() && search_ended_) {
        return "";
    }

    std::string res(solutions_queue_.front());
    solutions_queue_.pop_front();
    return res;
}

} /* namespace searchlight */

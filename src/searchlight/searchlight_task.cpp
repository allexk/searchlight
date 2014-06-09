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

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.task"));

void SearchlightTask::ResolveTask(const std::string &lib_name,
        const std::string &task_name) {
    // loading the task library
    void *lib_handle = dll_handler_.LoadDLL(lib_name);

    // look up the task
    task_ = (SLTaskFunc)dlsym(lib_handle, task_name.c_str());
    // We should check an error via dlerror, but NULL checking is fine
    if (!task_) {
        std::ostringstream err_msg;
        err_msg << "Cannot find an SL task function, name=" << task_name;
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << err_msg.str();
    }
}

std::string SearchlightTask::GetNextSolution() {
    boost::unique_lock<boost::mutex> lock(queue_mtx_);
    while (solutions_queue_.empty() && !search_ended_) {
        queue_cond_.wait(lock);
    }

    // Terminated on error? Throw it in the main thread.
    if (sl_error_) {
        LOG4CXX_INFO(logger, "SearchlightTask terminated on error!");
        sl_error_->raise();
    }

    if (solutions_queue_.empty() && search_ended_) {
        LOG4CXX_INFO(logger, "SearchlightTask finished");
        return "";
    }

    std::string res(solutions_queue_.front());
    solutions_queue_.pop_front();

    LOG4CXX_TRACE(logger, "A solution is found: " << res);

    return res;
}

const ConstChunk *SearchlightResultsArray::nextChunk(AttributeID attId,
        MemChunk& chunk) {
    // initialize the search if needed
    if (!sl_thread_) {
        LOG4CXX_INFO(logger, "Starting Searchlight Thread");
        sl_thread_ = new boost::thread(boost::ref(*sl_task_));
    }

    // block until next solution
    std::string sol = sl_task_->GetNextSolution();
    if (sol.empty()) {
        return NULL; //no more results
    } else {
        Coordinates pos(1, res_count_++);
        Address addr(attId, pos);
        chunk.initialize(this, &desc_, addr, 0 /* no compression */);

        // set iterator
        boost::shared_ptr<Query> query(Query::getValidQueryPtr(_query));
        boost::shared_ptr<ChunkIterator> dst = chunk.getIterator(query, 0);
        dst->setPosition(pos);

        // write value
        Value v;
        v.setString(sol.c_str());
        dst->writeItem(v);
        dst->flush();

        return &chunk;
    }
}

ArrayDesc SearchlightResultsArray::GetArrayDesc() {
    std::vector<AttributeDesc> attrs(1,
         AttributeDesc(AttributeID(0), "assignment", scidb::TID_STRING, 0, 0));
    std::vector<DimensionDesc> dims(1,
            DimensionDesc("num", 0, scidb::MAX_COORDINATE, 1, 0));
    return ArrayDesc("sl_result", attrs, dims);
}

} /* namespace searchlight */

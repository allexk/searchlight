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
#include "searchlight_messages.pb.h"
#include "validator.h"

#include <dlfcn.h>

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.task"));

void SearchlightTask::operator()() {
    try {
        bool work_expected = true;
        while (work_expected) {
            // Wait until something interesting happens
            Searchlight::Status sl_status;
            {
                std::unique_lock<std::mutex> lock(mtx_);
                sl_status = searchlight_.GetStatus();
                while (sl_status == Searchlight::Status::VOID) {
                    sl_cond_.wait(lock);
                    sl_status = searchlight_.GetStatus();
                }
            }

            // Either solve or exit
            switch (sl_status) {
                case Searchlight::Status::PREPARED:
                    searchlight_.Solve();
                    sl_status = searchlight_.GetStatus();
                    if (sl_status == Searchlight::Status::TERMINATED) {
                        LOG4CXX_DEBUG(logger, "Solver was terminated!");
                        work_expected = false;
                    }
                    break;
                case Searchlight::Status::FIN_SEARCH:
                case Searchlight::Status::FIN_VALID:
                case Searchlight::Status::TERMINATED:
                case Searchlight::Status::COMMITTED:
                    work_expected = false;
                    break;
                default:
                    assert(false);
                    throw SYSTEM_EXCEPTION(
                            SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                            << "Unexpected status in the SL loop!";
                    break;
            }
        }
    } catch (const scidb::Exception &ex) {
        /*
         * std::thread would call std::terminate(), so
         *  we have to take care of proper error reporting
         */
        sl_error_ = ex.copy();
        searchlight_.Terminate();
        queue_cond_.notify_one();
    } catch (const std::exception &e) {
        // Catch other C++ and library exceptions and translate them
        sl_error_ = (SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                SCIDB_LE_ILLEGAL_OPERATION) << e.what()).copy();
        searchlight_.Terminate();
        queue_cond_.notify_one();
    } catch (...) {
        sl_error_ = (SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                SCIDB_LE_ILLEGAL_OPERATION) <<
                "Unknown exception in SL!").copy();
        searchlight_.Terminate();
        queue_cond_.notify_one();
    }
}

void SearchlightTask::ReportIdleSolver() {
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    if (query->getCoordinatorID() == scidb::COORDINATOR_INSTANCE) {
        // We are at the coordinator -- handle the solver
        HandleIdleSolver(my_instance_id_);
    } else {
        // We are at a common instance -- send the message to the coordinator
        SearchlightMessenger::getInstance()->ReportIdleSolver(query);
    }
}

void SearchlightTask::ReportFinValidator() {
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    if (query->getCoordinatorID() == scidb::COORDINATOR_INSTANCE) {
        // We are at the coordinator -- handle the solver
        HandleFinValidator(my_instance_id_);
    } else {
        // We are at a common instance -- send the message to the coordinator
        SearchlightMessenger::getInstance()->ReportFinValidator(query);
    }
}

InstanceID SearchlightTask::GetHelpee() {
    /*
     * We have a queue of instances needing help. So, we just take the first
     * one. Then, we reinsert it in the back, so that we are able to send more
     * helpers there later.
     */
    InstanceID res = scidb::INVALID_INSTANCE;
    if (!distr_search_info_->helpees_.help_reqs_.empty()) {
        res = distr_search_info_->helpees_.help_reqs_.front();
        distr_search_info_->helpees_.help_reqs_.pop_front();
        distr_search_info_->helpees_.Add(res);
    }
    return res;
}

void SearchlightTask::HandleIdleSolver(InstanceID id) {
    std::unique_lock<std::mutex> lock(mtx_);
    assert(distr_search_info_->busy_solvers_.count(id));

    // The solver is not busy anymore and doesn't need help
    distr_search_info_->busy_solvers_.erase(id);
    distr_search_info_->idle_solvers_.insert(id);
    distr_search_info_->helpees_.Erase(id);

    // Check if the search is completely finished (not validation, though)
    if (distr_search_info_->busy_solvers_.empty()) {
        BroadcastFinishSearch();
    } else {
        CheckForHelpees();
    }
}

void SearchlightTask::CheckForHelpees() {
    while (!distr_search_info_->idle_solvers_.empty()) {
        const InstanceID helpee = GetHelpee();
        if (helpee != scidb::INVALID_INSTANCE) {
            auto helper_iter = distr_search_info_->idle_solvers_.begin();
            const InstanceID helper = *helper_iter;
            distr_search_info_->idle_solvers_.erase(helper_iter);
            distr_search_info_->busy_solvers_.insert(helper);
            DispatchHelper(helper, helpee);
        } else {
            break;
        }
    }
}

void SearchlightTask::DispatchHelper(InstanceID helper, InstanceID dest) {
    LOG4CXX_DEBUG(logger, "Dispatching help, dest=" << dest
            << ", helper=" << helper);
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    if (query->getCoordinatorID() == dest) {
        HandleHelper(helper);
    } else {
        SearchlightMessenger::getInstance()->DispatchHelper(
                query, helper, dest);
    }
}

void SearchlightTask::HandleHelper(InstanceID helper) {
    searchlight_.HandleHelper(helper);
}

void SearchlightTask::HandleAcceptHelper(InstanceID helper) {
    std::lock_guard<std::mutex> lock(mtx_);
    distr_search_info_->helpees_.Add(helper);
    CheckForHelpees();
}

void SearchlightTask::HandleEndOfSearch() {
    searchlight_.HandleEndOfSearch();
    sl_cond_.notify_one();
}

void SearchlightTask::HandleCommit() {
    searchlight_.HandleCommit();
    queue_cond_.notify_one();
}

void SearchlightTask::HandleFinValidator(InstanceID id) {
    std::unique_lock<std::mutex> lock(mtx_);
    assert(distr_search_info_->busy_solvers_.empty());
    assert(distr_search_info_->busy_validators_count_ > 0);

    // count the validator out
    distr_search_info_->busy_validators_count_--;

    if (distr_search_info_->busy_validators_count_ == 0) {
        BroadcastCommit();
    }
}

void SearchlightTask::BroadcastFinishSearch() {
    LOG4CXX_INFO(logger, "Broadcasting end-of-search...");
    // Local solver
    HandleEndOfSearch();
    // Remote solvers
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    SearchlightMessenger::getInstance()->BroadcastFinishSearch(query);
}

void SearchlightTask::BroadcastCommit() {
    LOG4CXX_INFO(logger, "Broadcasting commit...");
    // Local SL
    HandleCommit();
    // Remote SLs
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    SearchlightMessenger::getInstance()->BroadcastCommit(query);
}

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

void SearchlightTask::HandleRemoteSolution(InstanceID inst,
        const google::protobuf::Message *msg) {

    const SearchlightSolution *sol_msg =
            dynamic_cast<const SearchlightSolution *>(msg);
    if (!sol_msg) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Invalid message at the SL Task remote result function!";
    }

    std::vector<int64_t> vals;
    SearchlightMessenger::UnpackAssignment(sol_msg->solution(), vals, nullptr);
    ReportSolution(vals);
}

// Handles control messages
void SearchlightTask::HandleControlMessage(InstanceID inst,
        const google::protobuf::Message *msg) {
    // Get the message
    const SearchlightControl *control_msg =
            dynamic_cast<const SearchlightControl *>(msg);
    if (!control_msg) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Invalid message at the SL Task control handler function!";
    }

    switch (control_msg->type()) {
        case SearchlightControl::SEARCH_IDLE:
            LOG4CXX_DEBUG(logger, "Search is idle: id=" << inst);
            HandleIdleSolver(inst);
            break;
        case SearchlightControl::VALIDATOR_LOCAL_FIN:
            LOG4CXX_DEBUG(logger, "Validator finished: id=" << inst);
            HandleFinValidator(inst);
            break;
        case SearchlightControl::END_SEARCH:
            LOG4CXX_DEBUG(logger, "End-of-search request");
            HandleEndOfSearch();
            break;
        case SearchlightControl::COMMIT:
            LOG4CXX_DEBUG(logger, "Commit request");
            HandleCommit();
            break;
    }
}

// Handles balancing messages
void SearchlightTask::HandleBalanceMessage(InstanceID inst,
        const google::protobuf::Message *msg) {
    // Get the message
    const SearchlightBalance *balance_msg =
            dynamic_cast<const SearchlightBalance *>(msg);
    if (!balance_msg) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Invalid message at the SL Task balance handler function!";
    }

    switch (balance_msg->type()) {
        case SearchlightBalance::HELP_LOAD:
            LOG4CXX_DEBUG(logger, "Got a remote load for the solver...");
            HandleRemoteLoad(*balance_msg);
            break;
        case SearchlightBalance::REJECT_HELPER_HARD:
        case SearchlightBalance::REJECT_HELPER_SOFT:
        {
            const int inst_count = balance_msg->instance_size();
            assert(inst_count);
            std::vector<InstanceID> instances(inst_count);
            for (int i = 0; i < inst_count; i++) {
                instances[i] = balance_msg->instance(i);
            }
            HandleRejectHelp(inst, instances,
                 balance_msg->type() == SearchlightBalance::REJECT_HELPER_HARD);
            break;
        }
        case SearchlightBalance::HELPER_DISPATCH:
            assert(balance_msg->instance_size());
            for (int i = 0; i < balance_msg->instance_size(); i++) {
                HandleHelper(balance_msg->instance(i));
            }
            break;
        case SearchlightBalance::CANDIDATE_FORWARD:
            LOG4CXX_DEBUG(logger, "Got forwards for the validator ...");
            HandleForwards(*balance_msg, inst);
            break;
        case SearchlightBalance::BALANCE_RESULT:
            searchlight_.GetValidator().HandleForwardResult(
                    balance_msg->id(), balance_msg->result());
            break;
        case SearchlightBalance::ACCEPT_HELP:
            assert(balance_msg->instance_size());
            for (int i = 0; i < balance_msg->instance_size(); i++) {
                HandleAcceptHelper(balance_msg->instance(i));
            }
            break;
    }
}

void SearchlightTask::HandleRemoteLoad(const SearchlightBalance &msg) {
    // Prepare the load
    LiteAssignmentVector load(msg.load_size());
    for (int i = 0; i < msg.load_size(); i++) {
        SearchlightMessenger::UnpackAssignment(msg.load(i), load[i]);
    }

    // Put the load into the solver
    std::lock_guard<std::mutex> lock(mtx_);
    searchlight_.PrepareHelper(load);
    sl_cond_.notify_one();
}

void SearchlightTask::HandleForwards(const SearchlightBalance &msg,
        InstanceID src) {
    // Prepare the load
    LiteAssignmentVector load(msg.load_size());
    for (int i = 0; i < msg.load_size(); i++) {
        SearchlightMessenger::UnpackAssignment(msg.load(i), load[i]);
    }

    searchlight_.GetValidator().AddRemoteCandidates(load, src, msg.id());
}


void SearchlightTask::HandleRejectHelp(InstanceID src,
        const std::vector<InstanceID> &helpers, bool hard) {
    // logging
    if (logger->isDebugEnabled()) {
        std::ostringstream deb_str;
        deb_str << "Help was rejected: helpee=" << src << "hard="
                << hard << ", helpers={";
        std::copy(helpers.begin(), helpers.end(),
                std::ostream_iterator<InstanceID>(deb_str, ", "));
        deb_str << "}";
        logger->debug(deb_str.str(), LOG4CXX_LOCATION);
    }

    // If it's a hard reject we need to disable help for src solver
    if (hard) {
        std::lock_guard<std::mutex> lock(mtx_);
        distr_search_info_->helpees_.Erase(src);
    }

    // Handle a rejected helper as idle -- there is no difference
    for (auto id: helpers) {
        HandleIdleSolver(id);
    }
}

void SearchlightTask::RejectHelp(const std::vector<InstanceID> &helpers,
        bool hard) {
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    if (query->getCoordinatorID() == scidb::COORDINATOR_INSTANCE) {
        HandleRejectHelp(my_instance_id_, helpers, hard);
    } else {
        // We are at a common instance -- send helpers back to the coordinator
        SearchlightMessenger::getInstance()->RejectHelp(query, helpers, hard);
    }
}

void SearchlightTask::ReportSolution(const std::vector<int64_t> &values) {
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    if (query->getCoordinatorID() == scidb::COORDINATOR_INSTANCE) {
        // We are at the coordinator -- add the solution to the queue
        const std::string sol = searchlight_.SolutionToString(values);
        std::lock_guard<std::mutex> lock(mtx_);
        solutions_queue_.push_back(sol);
        queue_cond_.notify_one();
    } else {
        // We are at a common instance -- send the solution to the coordinator
        SearchlightMessenger::getInstance()->SendSolution(query, values);
    }
}

std::string SearchlightTask::GetNextSolution() {
    std::unique_lock<std::mutex> lock(mtx_);
    while (solutions_queue_.empty() && ExpectingResults() && !sl_error_) {
        queue_cond_.wait(lock);
    }

    // Terminated on error? Throw it in the main thread.
    if (sl_error_) {
        LOG4CXX_INFO(logger, "SearchlightTask terminated on error!");
        sl_error_->raise();
    }

    if (solutions_queue_.empty()) {
        assert(!ExpectingResults());
        LOG4CXX_INFO(logger, "SearchlightTask finished");
        return "";
    }

    std::string res(solutions_queue_.front());
    solutions_queue_.pop_front();

    LOG4CXX_TRACE(logger, "A solution is found: " << res);

    return res;
}

void SearchlightTask::DispatchWork(const LiteAssignmentVector &work,
        InstanceID solver) {
    // For now a helper is always a remote instance.
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    SearchlightMessenger::getInstance()->DispatchWork(query, work, solver);

    // Next, we notify the coordinator about acceptance
    if (query->getCoordinatorID() == scidb::COORDINATOR_INSTANCE) {
        HandleAcceptHelper(solver);
    } else {
        SearchlightMessenger::getInstance()->AcceptHelp(query, solver);
    }
}

void SearchlightTask::ForwardCandidates(const LiteAssignmentVector &cands,
        InstanceID dest, int forw_id) const {
    // For now a helper is always a remote instance.
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    SearchlightMessenger::getInstance()->ForwardCandidates(query, cands, dest,
            forw_id);
}

void SearchlightTask::SendBalanceResult(InstanceID dest, int id,
        bool result) const {
    // For now a helper is always a remote instance.
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    SearchlightMessenger::getInstance()->SendBalanceResult(query, dest, id,
            result);
}

const ConstChunk *SearchlightResultsArray::nextChunk(AttributeID attId,
        MemChunk& chunk) {
    /*
     * We wait for results only at the coordinator, since it collects all of
     * them. Other instances will use SL Messenger to transfer local results,
     * so they should return EOF for the array immediately.
     */
    {
        boost::shared_ptr<Query> query(Query::getValidQueryPtr(_query));
        if (query->getCoordinatorID() != scidb::COORDINATOR_INSTANCE) {
            return NULL;
            // The searchlight and validator threads will continue working
        }
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

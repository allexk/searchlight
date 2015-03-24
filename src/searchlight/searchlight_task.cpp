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

/*
 * There is a compile warning here about an always-true comparison in
 * json_parser_write inside the header. Actually, the warning is
 * valid for char-based strings, but harmless. See boost bug #5598
 */
#if defined(__GNUC__) && (((__GNUC__ * 100) + __GNUC_MINOR__) >= 406)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wtype-limits"
#endif
#include <boost/property_tree/json_parser.hpp>
#if defined(__GNUC__) && (((__GNUC__ * 100) + __GNUC_MINOR__) >= 406)
#pragma GCC diagnostic pop
#endif

#include <boost/lexical_cast.hpp>

namespace {
/*
 * Read and parses a JSON array of T at the specified path.
 */
template<typename T>
std::vector<T> ReadJSONArray(
        searchlight::SearchlightConfig &sl_config,
        const std::string &path) {
    std::vector<T> res;
    try {
        for (const auto &node: sl_config.get_child(path)) {
            const T val = node.second.get_value<T>();
            res.push_back(val);
        }
    } catch (const boost::property_tree::ptree_bad_path &) {
        // ignore errors
    }
    return res;
}

/*
 * Read and parses a JSON map at the specified path.
 */
template<typename K, typename V>
std::map<K, V> ReadJSONMap(
        searchlight::SearchlightConfig &sl_config,
        const std::string &path) {
    std::map<K, V> res;
    try {
        for (const auto &node: sl_config.get_child(path)) {
            const K key = boost::lexical_cast<K>(node.first);
            const V val = boost::lexical_cast<V>(node.second.data());
            res[key] = val;
        }
    } catch (const boost::property_tree::ptree_bad_path &) {
        // ignore errors
    }
    return res;
}

/*
 * Solver id that is not equivalent to any valid id.
 */
const uint64_t INVALID_SOLVER_ID = ~uint64_t(0);

inline scidb::InstanceID GetInstanceFromSolverID(uint64_t id) {
    // Just ignore low 32 bits
    return id >> 32;
}

inline uint32_t GetLocalIDFromSolverID(uint64_t id) {
    // Just ignore high 32 bits
    return uint32_t(id);
}

inline uint64_t CreateSolverID(scidb::InstanceID inst, uint32_t ord) {
    return (uint64_t(inst) << 32) + ord;
}

} /* namespace <anonymous> */

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.task"));

// The logger for results
static log4cxx::LoggerPtr result_logger(
        log4cxx::Logger::getLogger("searchlight.result"));

SearchlightTask::SearchlightTask(const std::string &library_name,
        const std::string &task_name, const std::string &config_file_name,
        ArrayPtr &data, const ArrayPtrVector &samples,
        const boost::shared_ptr<Query> &query) :
            query_instance_count_(query->getInstancesCount()),
            my_instance_id_(query->getInstanceID()),
            searchlight_(*this, dll_handler_),
            query_(query) {

    // Read the config
    ReadConfig(config_file_name);

    // Fill in distributed search info
    if (query->getCoordinatorID() == scidb::COORDINATOR_INSTANCE) {
        distr_search_info_.reset(new DistributedSearchInfo);
        for (size_t i = 0; i < active_solver_instances_.size(); i++) {
            // Each instance might contain several solvers
            const int solv_count = active_solver_num_[i];
            uint64_t solv_id = CreateSolverID(active_solver_instances_[i], 0);
            for (int s = 0; s < solv_count; s++, solv_id++) {
                distr_search_info_->helpees_.Add(solv_id);
                distr_search_info_->busy_solvers_.insert(solv_id);
            }
        }

        /*
         *  We need the list of active instances here. An instance is active if
         *  it contains an active solver or validator. What's important is that
         *  these instances will contain a validator (active or just
         *  a forwarder), and we need to count all of them.
         *
         *  We use set here for de-duplication.
         */
        std::unordered_set<InstanceID> active_instances{
            active_solver_instances_.begin(),
            active_solver_instances_.end()};
        active_instances.insert(active_validators_.begin(),
                active_validators_.end());
        distr_search_info_->busy_validators_count_ = active_instances.size();
    }

    ResolveTask(library_name, task_name);
    searchlight_.RegisterArray(data, samples);

    // Are we using dynamic scheduling for helpers/solvers?
    dynamic_scheduling_ =
            config_.get("searchlight.dynamic_scheduling", 1);

    // Set distribution update frequence for the data array
    const int distr_update_freq =
            config_.get("balance.map_update_frequency", 5);
    SearchlightMessenger::getInstance()->SetDistributedMapUpdateFrequency(
            query, data->getArrayDesc().getName(), distr_update_freq);

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

    // Prepare Searchlight
    if (InstanceActive(my_instance_id_)) {
        int solver_count = GetSolverNum(my_instance_id_);
        if (solver_count == 0) {
            // Create a dummy solver to hold the model (needed by validator)
            solver_count = 1;
        }
        searchlight_.Prepare(task_name, task_,
                CreateSolverID(my_instance_id_, 0), solver_count);
    }
}

size_t SearchlightTask::GetGlobalOrdinalSolverId(uint64_t id) const {
    // Instance ordinal: high 32 bits
    const InstanceID inst = GetInstanceFromSolverID(id);
    const auto iter = std::find(active_solver_instances_.begin(),
            active_solver_instances_.end(), inst);
    assert(iter != active_solver_instances_.end());
    const int inst_ord = iter - active_solver_instances_.begin();

    // Global solver ordinal
    size_t res = 0;
    for (int i = 0; i < inst_ord; i++) {
        res += active_solver_num_[i];
    }
    res += GetLocalIDFromSolverID(id);

    return res;
}

void SearchlightTask::ReadConfig(const std::string &file_name) {
    try {
        boost::property_tree::read_json(file_name, config_);
    } catch (const boost::property_tree::json_parser_error &e) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << e.what();
    }

    if (logger->isDebugEnabled()) {
        std::ostringstream deb;
        deb << "Config for the task follows:\n";
        boost::property_tree::write_json(deb, config_, true);
        logger->debug(deb.str(), LOG4CXX_LOCATION);
    }

    // Active solvers
    auto solvers_setup = ReadJSONMap<InstanceID, int>(config_,
            "setup.solvers");
    if (solvers_setup.empty()) {
        for (int i = 0; i < query_instance_count_; i++) {
            solvers_setup[i] = 1;
        }
    }
    active_solver_instances_.reserve(solvers_setup.size());
    active_solver_num_.reserve(solvers_setup.size());
    for (auto it = solvers_setup.begin(); it != solvers_setup.end(); ++it) {
        if (it->second > 0) {
            active_solver_instances_.push_back(it->first);
            active_solver_num_.push_back(it->second);
        }
    }

    // Active validators
    active_validators_ = ReadJSONArray<InstanceID>(config_, "setup.validators");
    if (active_validators_.empty()) {
        active_validators_.resize(query_instance_count_);
        for (int i = 0; i < query_instance_count_; i++) {
            active_validators_[i] = i;
        }
    }
}

void SearchlightTask::ReportIdleSolver(uint64_t solver_id, bool deferable) {
    if (dynamic_scheduling_ && deferable) {
        /*
         * For now we prefer validator helpers to pending loads. So, we first
         * check if the validator wants to run a helper, even just for a
         * round-trip to the coordinator to report idleness and get a new
         * load.
         *
         * Another way might be to perform: FreeThread(); ReserveThread() to
         * give pending loads the ability to cease the thread.
         */
        bool persistent_helper;
        const bool started = searchlight_.GetValidator().AddValidatorHelper(
                &persistent_helper);
        if (started && persistent_helper) {
            // If the helper is persistent, don't report the solver as idle yet
            std::lock_guard<std::mutex> lock{mtx_};
            pending_idle_solvers_.insert(solver_id);
            return;
        }
    }
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    if (query->getCoordinatorID() == scidb::COORDINATOR_INSTANCE) {
        // We are at the coordinator -- handle the solver
        HandleIdleSolver(solver_id);
    } else {
        // We are at a common instance -- send the message to the coordinator
        SearchlightMessenger::getInstance()->ReportIdleSolver(query, solver_id);
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

uint64_t SearchlightTask::GetHelpee() {
    /*
     * We have a queue of instances needing help. So, we just take the first
     * one. Then, we reinsert it in the back, so that we are able to send more
     * helpers there later.
     */
    uint64_t res = INVALID_SOLVER_ID;
    if (!distr_search_info_->helpees_.help_reqs_.empty()) {
        res = distr_search_info_->helpees_.help_reqs_.front();
        distr_search_info_->helpees_.Erase(res);
        distr_search_info_->helpees_.Add(res);
    }
    return res;
}

void SearchlightTask::HandleIdleSolver(uint64_t id) {
    std::unique_lock<std::mutex> lock(mtx_);
    assert(distr_search_info_->busy_solvers_.count(id));

    // The solver is not busy anymore and doesn't need help
    distr_search_info_->busy_solvers_.erase(id);
    distr_search_info_->idle_solvers_.insert(id);
    distr_search_info_->helpees_.Erase(id);

    /*
     * Check if the search is completely finished. Validator might be
     * still going.
     *
     * Here, we try to unlock the mutex ASAP. This is not only a performance
     * consideration. There might be a call path: sl -> sl_task -> sl that
     * leads to double locking. Since such paths are rare, we don't use
     * a recursive mutex due to the cost.
     */
    const bool search_finished = distr_search_info_->busy_solvers_.empty();
    lock.unlock();

    if (search_finished) {
        BroadcastFinishSearch();
    } else {
        CheckForHelpees();
    }
}

void SearchlightTask::CheckForHelpees() {
    std::unique_lock<std::mutex> lock(mtx_);
    while (!distr_search_info_->idle_solvers_.empty()) {
        const InstanceID helpee = GetHelpee();
        if (helpee != INVALID_SOLVER_ID) {
            auto helper_iter = distr_search_info_->idle_solvers_.begin();
            const InstanceID helper = *helper_iter;
            distr_search_info_->idle_solvers_.erase(helper_iter);
            distr_search_info_->busy_solvers_.insert(helper);

            /*
             * Unlock here, since the call might return to the task via SL.
             * For example, the helper goes to the same instance and is
             * immediately rejected. Reject requires another mutex lock
             * in the same thread.
             *
             * See also the comment in HandleIdleSolver().
             */
            lock.unlock();
            DispatchHelper(helper, helpee);
            lock.lock();
        } else {
            break;
        }
    }
}

void SearchlightTask::DispatchHelper(uint64_t helper, uint64_t dest) {
    LOG4CXX_DEBUG(logger, "Dispatching help, dest=0x" << std::hex << dest
            << ", helper=0x" << helper);
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    const InstanceID dest_inst = GetInstanceFromSolverID(dest);
    if (my_instance_id_ == dest_inst) {
        HandleHelper(helper, dest);
    } else {
        SearchlightMessenger::getInstance()->DispatchHelper(
                query, helper, dest, dest_inst);
    }
}

void SearchlightTask::HandleHelper(uint64_t helper, uint64_t helpee) {
    /*
     * Special conditions (like "pending idle") are only relevant for VOID
     * solvers. For others, we let them handle the helper themselves.
     */
    if (dynamic_scheduling_ &&
            searchlight_.GetSolverStatus(GetLocalIDFromSolverID(helpee)) ==
                    SearchlightSolver::Status::VOID) {
        std::unique_lock<std::mutex> lock{mtx_};
        /*
         * First, we want to see if the solver has a pending assignment. If so,
         * we might try to off-load it to the helper.
         */
        const auto iter = pending_assgns_.find(helpee);
        if (iter != pending_assgns_.end()) {
            /*
             * We have a pending assignment -- give it away to the helper
             * TODO: Think about time-stamps to avoid unnecessary load movements.
             */
            if (GetInstanceFromSolverID(helper) ==
                    GetInstanceFromSolverID(helpee)) {
                // The same instance -- no gain, since we don't have threads
                lock.unlock();
                RejectHelp({helper}, helpee, false);
                LOG4CXX_DEBUG(logger, "Got the same instance helper for a"
                        "pending assignment, reject...");
            } else {
                LiteAssignmentVector load_copy{std::move(iter->second)};
                pending_assgns_.erase(iter);
                lock.unlock();
                DispatchWork(load_copy, helper);
                ReportIdleSolver(helpee, false);
                LOG4CXX_DEBUG(logger, "Got a helper for a pending assignment, "
                        "forwarding the load...");
            }
            return;
        }

        /*
         * Secondly, we check if the solver is pending idle. In this case it's
         * still busy at the coordinator. We'll make a hard help reject then,
         * to make sure no helpers arrive until the solver reports itself idle.
         */
        if (pending_idle_solvers_.count(helpee)) {
            lock.unlock();
            RejectHelp({helper}, helpee, true);
            LOG4CXX_DEBUG(logger, "Got a helper for a pending idle solver, "
                    "hard reject, id=0x" << std::hex << helpee);
            return;
        }
    }
    searchlight_.HandleHelper(helper, GetLocalIDFromSolverID(helpee));
}

void SearchlightTask::HandleAcceptHelper(uint64_t helper) {
    std::unique_lock<std::mutex> lock(mtx_);
    // Check if the message came too late; solver could report idleness already
    if (distr_search_info_->busy_solvers_.find(helper) !=
            distr_search_info_->busy_solvers_.end()) {
        distr_search_info_->helpees_.Add(helper);
        lock.unlock();
        CheckForHelpees();
    }
}

void SearchlightTask::HandleEndOfSearch() {
    if (InstanceActive(my_instance_id_)) {
        searchlight_.HandleEndOfSearch();
    }
}

void SearchlightTask::HandleCommit() {
    if (InstanceActive(my_instance_id_)) {
        searchlight_.HandleCommit();
    }
    queue_cond_.notify_one();
}

void SearchlightTask::HandleFinValidator(InstanceID id) {
    std::unique_lock<std::mutex> lock(mtx_);
    assert(distr_search_info_->busy_solvers_.empty());
    assert(distr_search_info_->busy_validators_count_ > 0);

    // count the validator out
    distr_search_info_->busy_validators_count_--;

    if (distr_search_info_->busy_validators_count_ == 0) {
        lock.unlock();
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
            LOG4CXX_DEBUG(logger, "Search is idle: id=0x" << std::hex
                    << control_msg->id(0));
            HandleIdleSolver(control_msg->id(0));
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
            assert(balance_msg->id_size());
            LOG4CXX_DEBUG(logger, "Got a remote load for the solver, id=0x"
                    << std::hex << balance_msg->id(0));
            HandleRemoteLoad(*balance_msg, balance_msg->id(0));
            break;
        case SearchlightBalance::REJECT_HELPER_HARD:
        case SearchlightBalance::REJECT_HELPER_SOFT:
        {
            // First id is the sender-solver
            const int solv_count = balance_msg->id_size() - 1;
            assert(solv_count);
            std::vector<uint64_t> ids(solv_count);
            for (int i = 0; i < solv_count; i++) {
                ids[i] = balance_msg->id(i + 1);
            }
            HandleRejectHelp(balance_msg->id(0), ids,
                 balance_msg->type() == SearchlightBalance::REJECT_HELPER_HARD);
            break;
        }
        case SearchlightBalance::HELPER_DISPATCH:
            assert(balance_msg->id_size() >= 2);
            // First id is helpee, the rest are helpers
            for (int i = 0; i < balance_msg->id_size() - 1; i++) {
                HandleHelper(balance_msg->id(i + 1), balance_msg->id(0));
            }
            break;
        case SearchlightBalance::CANDIDATE_FORWARD:
            LOG4CXX_TRACE(logger,
                    "Got forwards for the validator from inst=" << inst);
            HandleForwards(*balance_msg, inst);
            break;
        case SearchlightBalance::BALANCE_RESULT:
            searchlight_.GetValidator().HandleForwardResult(
                    balance_msg->id(0), balance_msg->result());
            break;
        case SearchlightBalance::ACCEPT_HELP:
            assert(balance_msg->id_size());
            for (int i = 0; i < balance_msg->id_size(); i++) {
                LOG4CXX_DEBUG(logger, "Helper was accepted: helper=0x"
                        << std::hex << balance_msg->id(i));
                HandleAcceptHelper(balance_msg->id(i));
            }
            break;
        case SearchlightBalance::VALIDATOR_INFO:
            assert(balance_msg->id_size() && balance_msg->id(0) >= 0);
            assert(std::find(active_validators_.begin(),
                    active_validators_.end(), inst) !=
                            active_validators_.end());

            const size_t cands_num = balance_msg->id(0);
            const size_t val_ord = std::distance(active_validators_.begin(),
                    std::find(active_validators_.begin(),
                            active_validators_.end(), inst));
            if (InstanceActive(my_instance_id_)) {
                LOG4CXX_DEBUG(logger, "Validator broadcasted non-simulated "
                        "candidates number: val=" << val_ord <<
                        ", cands=" << cands_num);
                searchlight_.GetValidator().UpdateCandsInfo(cands_num, val_ord);
            }
            break;
    }
}

void SearchlightTask::HandleRemoteLoad(const SearchlightBalance &msg,
        uint64_t helper) {
    // Prepare the load
    LiteAssignmentVector load(msg.load_size());
    for (int i = 0; i < msg.load_size(); i++) {
        SearchlightMessenger::UnpackAssignment(msg.load(i), load[i]);
    }

    // Put or store the load
    DispatchOrStoreLoad(load, helper);
}

void SearchlightTask::HandleForwards(const SearchlightBalance &msg,
        InstanceID src) {
    // Prepare the load
    LiteAssignmentVector load(msg.load_size());
    Int64Vector zones(msg.load_size()); // 0 is the default zone
    for (int i = 0; i < msg.load_size(); i++) {
        const auto &load_i = msg.load(i);
        SearchlightMessenger::UnpackAssignment(load_i, load[i]);
        if (load_i.aux_info_size()) {
            zones[i] = load_i.aux_info(0);
        }
    }
    searchlight_.GetValidator().AddRemoteCandidates(load, zones, src,
            msg.id(0));
}

void SearchlightTask::HandleRejectHelp(uint64_t src,
        const std::vector<uint64_t> &helpers, bool hard) {
    // logging
    if (logger->isDebugEnabled()) {
        std::ostringstream deb_str;
        deb_str << "Help was rejected: helpee=0x" << std::hex << src
                << ", hard=" << hard << ", helpers={0x";
        std::copy(helpers.begin(), helpers.end(),
                std::ostream_iterator<InstanceID>(deb_str, ", 0x"));
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

void SearchlightTask::RejectHelp(const std::vector<uint64_t> &helpers,
        uint64_t solver_id, bool hard) {
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    if (query->getCoordinatorID() == scidb::COORDINATOR_INSTANCE) {
        HandleRejectHelp(solver_id, helpers, hard);
    } else {
        // We are at a common instance -- send helpers back to the coordinator
        SearchlightMessenger::getInstance()->RejectHelp(query, helpers,
                solver_id, hard);
    }
}

void SearchlightTask::ReportSolution(const std::vector<int64_t> &values) {
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    if (query->getCoordinatorID() == scidb::COORDINATOR_INSTANCE) {
        // We are at the coordinator -- add the solution to the queue
        const std::string sol = searchlight_.SolutionToString(values);
        LOG4CXX_INFO(result_logger, sol);
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
        queue_cond_.wait_for(lock, std::chrono::seconds(10));
        /*
         * We wait for 10 seconds to check errors one in a while. Such
         * errors might occur when a worker sends mtError, but the
         * coordinator doesn't do anything except sitting here and
         * answering messages.
         */
        Query::getValidQueryPtr(query_);
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

bool SearchlightTask::PendingSolverJobs(bool check_idle_solvers) {
    if (!dynamic_scheduling_) {
        return false;
    }

    std::unique_lock<std::mutex> lock{mtx_};
    const bool res = !pending_assgns_.empty();
    if (check_idle_solvers && !pending_idle_solvers_.empty()) {
        const auto iter = pending_idle_solvers_.begin();
        const uint64_t idle_solver = *iter;
        pending_idle_solvers_.erase(iter);
        LOG4CXX_DEBUG(logger, "Releasing pending idle solver, id=0x" <<
                std::hex << idle_solver);
        lock.unlock();
        ReportIdleSolver(idle_solver, false);
    }
    return res;
}

void SearchlightTask::FreeThread() {
    if (dynamic_scheduling_) {
        std::unique_lock<std::mutex> lock{mtx_};
        if (!pending_assgns_.empty()) {
            // Have a pending load
            const auto load_it = pending_assgns_.begin();
            LiteAssignmentVector load_copy{std::move(load_it->second)};
            const uint64_t helper = load_it->first;
            pending_assgns_.erase(load_it);
            lock.unlock();

            searchlight_.PrepareHelper(load_copy,
                    GetLocalIDFromSolverID(helper));
            LOG4CXX_INFO(logger, "Dispatched a pending assignment...");
        } else {
            ++threads_available_;
        }
    }
}

void SearchlightTask::DispatchWork(const LiteAssignmentVector &work,
        uint64_t dest_solver) {
    // First, dispatch work
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    const InstanceID dest_inst = GetInstanceFromSolverID(dest_solver);
    if (my_instance_id_ == dest_inst) {
        // helper is a local solver
        LiteAssignmentVector work_copy{work};
        DispatchOrStoreLoad(work_copy, dest_solver);
    } else {
        SearchlightMessenger::getInstance()->DispatchWork(query, work,
                dest_solver, dest_inst);
    }

    // Next, we notify the coordinator about acceptance
    if (query->getCoordinatorID() == scidb::COORDINATOR_INSTANCE) {
        HandleAcceptHelper(dest_solver);
    } else {
        SearchlightMessenger::getInstance()->AcceptHelp(query, dest_solver);
    }
}

void SearchlightTask::DispatchOrStoreLoad(LiteAssignmentVector &load,
        uint64_t dest_solver) {
    if (dynamic_scheduling_) {
        std::lock_guard<std::mutex> lock{mtx_};
        if (threads_available_ > 0) {
            // Reserve a thread
            threads_available_--;
        } else {
            // No threads -- have to save the load for later
            assert(pending_assgns_.count(dest_solver) == 0);
            pending_assgns_.emplace(dest_solver, std::move(load));
            LOG4CXX_INFO(logger,
                    "No threads available: saving the remote load...");
            return;
        }
    }

    // Put the load into the solver
    searchlight_.PrepareHelper(load, GetLocalIDFromSolverID(dest_solver));
}

void SearchlightTask::ForwardCandidates(const LiteAssignmentVector &cands,
        const std::vector<int64_t> &zones,
        InstanceID dest, uint64_t forw_id) const {
    // For now a helper is always a remote instance.
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    SearchlightMessenger::getInstance()->ForwardCandidates(query, cands,
            zones, dest, forw_id);
}

void SearchlightTask::SendBalanceResult(InstanceID dest, uint64_t id,
        bool result) const {
    // For now a helper is always a remote instance.
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    SearchlightMessenger::getInstance()->SendBalanceResult(query, dest, id,
            result);
}

void SearchlightTask::BroadcastValidatorInfo(size_t cands_num) const {
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    SearchlightMessenger::getInstance()->BroadcastValidatorInfo(query,
            cands_num);
}

void SearchlightTask::HandleSearchlightError(
        const scidb::Exception::Pointer &error) {
    std::lock_guard<std::mutex> lock{mtx_};
    if (!sl_error_) {
        sl_error_ = error;
        Terminate();
        try {
            // We should try to notify the query
            const auto query = Query::getValidQueryPtr(query_);
            query->handleError(error);
            if (query->getCoordinatorID() != scidb::COORDINATOR_INSTANCE) {
                // Send to the coordinator
                auto error_msg = scidb::makeErrorMessageFromException(
                        *sl_error_, query->getQueryID());
                NetworkManager::getInstance()->sendMessage(
                        query->getPhysicalCoordinatorID(), error_msg);
                LOG4CXX_INFO(logger, "Notified coordinator about the error");
            }
        } catch (const scidb::Exception &) {
            // Query might have already been deallocated
            LOG4CXX_INFO(logger, "Reported error, but query is already "
                    "deallocating...");
        }
        queue_cond_.notify_one();
    }
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

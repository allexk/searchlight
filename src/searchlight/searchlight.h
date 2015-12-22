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
 * @file searchlight.h
 *
 * This is the main entry point for the search process. It contains a number of
 * tools to guide and monitor search, and also provides a number of other
 * useful tools, like estimation API.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_SEARCHLIGHT_H_
#define SEARCHLIGHT_SEARCHLIGHT_H_

#include "ortools_inc.h"
#include "array_desc.h"

#include <system/Config.h>

#include <dlfcn.h>
#include <boost/thread.hpp>
#include <boost/make_shared.hpp>
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>

namespace searchlight {

class Validator;
class SearchlightCollector;
class SearchlightTask;
class Searchlight;

/**
 * Type of the task function (called from the library).
 *
 * The second parameter is the id of the solver.
 */
typedef void (*SLTaskFunc)(Searchlight *, uint32_t);

/**
 * This class manages supplementary DLL resources required by Searchlight.
 * It is guaranteed to close all opened libraries at the end, on destruction.
 * It also serves as a map of names->dll handles, thus allowing retrieval of
 * DLLs via names, avoiding duplicate open calls.
 *
 * The need of such a manager is also emphasized by the necessity of closing
 * the libraries at the very end. For example, if you create an object with
 * the definition in a DLL and later try to delete it after closing the
 * library, you will have trouble with calling the destructor. That is why
 * the handler should be deleted at the very end.
 */
class DLLHandler : private boost::noncopyable {
public:
    /**
     * Constructs a new DLL handler. It assumes the default DLL directory
     * to be the SciDb plugin's directory.
     */
    DLLHandler() :
        dlls_dir_(scidb::Config::getInstance()->
                getOption<std::string>(scidb::CONFIG_PLUGINSDIR)) {}

    /**
     * Destructor. Closes all opened DLLs.
     */
    ~DLLHandler() {
        for (auto &name_lib: dlls_) {
            dlclose(name_lib.second);
        }

    }

    /**
     * Loads the specified DLL into memory. If the library was loaded before,
     * it returns the same handle. If the loading is impossible for some
     * reason, it throws a SciDb system exception.
     *
     * @param name the name of the library, without the suffix or prefix
     * @return the DLL handle (never nullptr)
     */
    void *LoadDLL(const std::string &name) {
        auto it = dlls_.find(name);
        if (it == dlls_.end()) {
            //std::string full_name = dlls_dir_ + "/lib" + name + ".so";
            std::string full_name = std::string{"lib"} + name + ".so";
            void *dll_handle= dlopen(full_name.c_str(), RTLD_LAZY | RTLD_LOCAL);
            if (!dll_handle) {
                std::ostringstream err_msg;
                err_msg << "Cannot load the task library: name=" <<
                        full_name << ", reason=" << dlerror();
                throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                        SCIDB_LE_ILLEGAL_OPERATION) << err_msg.str();
            }

            dlls_[name] = dll_handle;
            return dll_handle;
        } else {
            return it->second;
        }
    }

private:
    // Default DLL directory
    std::string dlls_dir_;

    // Map: DLL name -> DLL handle
    std::map<std::string, void *> dlls_;
};

/**
 * This is a monitor that catches complete, but approximate, solutions
 * (assignments) and passed them along to the Validator.
 */
class ValidatorMonitor : public SolutionCollector {
public:
    /**
     * Creates a new validator monitor. This monitor looks for complete
     * assignments during the search and passed them along to the validator.
     *
     * @param vars a vector of decision variables (externally managed)
     * @param solver the main solver
     * @param validator validator to connect
     */
    ValidatorMonitor(const IntVarVector &vars, Solver *solver,
            Validator &validator) :
        SolutionCollector(solver),
        validator_(validator),
        vars_(vars) {
        Add(vars);
    }

    /**
     * This function is called at a leaf of the search tree. At this point
     * a leaf is accepted as being a solution. This function checks if
     * it is a complete assignment and passed it along to the validator.
     *
     * @return true if we want to continue after the leaf; false otherwise
     */
    virtual bool AtSolution();

    /**
     * Returns the number of candidates encountered. This monitor tracks them,
     * since the number of leaves is not transfered between searches in the
     * or-tools solver. Thus, it creates problems when nested searches are
     * used.
     *
     * @return the number of candidates
     */
    int64_t CandidatesNumber() const {
        return candidates_;
    }

private:
    // The validator to pass the solution to
    Validator &validator_;

    // The vector of vars (managed outside)
    const IntVarVector vars_;

    // Candidates encountered
    int64_t candidates_ = 0;
};

/**
 * The type for a UDF function creator. It produces an or-tools IntExpr
 * representing the function. Takes as parameters: the solver to use
 * with, the adapter for accessing data, a vector of variables to work
 * with and integer parameters.
 */
typedef IntExpr *(* UDFFunctionCreator)(Solver *, AdapterPtr,
        const std::vector<IntVar *> &, const std::vector<int64> &);

/*
 * This class encapsulates a single solver thread-instance. The model,
 * monitors, heuristics go here. The Searchlight class provides higher
 * level first-second level interface.
 */
class SearchlightSolver : private boost::noncopyable {
public:
    /**
     * Status of Searchlight solver.
     */
    enum class Status {
        VOID,      //!< VOID Non-initialized state
        PREPARED,  //!< PREPARED Main solver prepared
        SEARCHING, //!< SEARCHING Main solver is searching (suring Solve())
        FINISHED,  //!< FINISHED Main solver finished completely
        TERMINATED //!< TERMINATED Abnormal termination on error
    };

    /**
     * Constructs a new solver.
     *
     * @param sl searchlight instance
     * @param id solver's id (instance + global ordinal)
     * @param sl_name name of the model (search)
     */
    SearchlightSolver(Searchlight &sl, uint64_t id,
            const std::string &sl_name) :
        sl_(sl),
        solver_(sl_name + "_" + std::to_string(id)),
        db_(nullptr),
        vars_leaf_(&solver_),
        status_(Status::VOID),
        id_(id) {}

    /**
     * Destructor.
     */
    ~SearchlightSolver();

    /**
     * The main solve method that starts the search.
     *
     * It is assumed that the
     * solver was prepared by calling Prepare() function.
     */
    void Solve();

    /**
     * Set the task for execution. This method does not start the
     * search. The search itself is started by the engine.
     *
     * This method takes "primary" and "secondary" decision variables.
     * "Primary" variables basically define the location of the object of
     * search (e.g., the left-most corner of a rectangle). "Secondary"
     * variables convey additional information about the object, like perimeter,
     * length and so on. Essentially, primary variables are used to determine
     * search distribution across multiple instances.
     *
     * @param primary_vars "primary" decision variables
     * @param secondary_vars "secondary" decision variables
     * @param db the decision builder (search heuristic)
     * @param vars the decision variables
     * @param monitors monitors, if required (can be empty)
     */
    void SetTask(const IntVarVector &primary_vars,
            const IntVarVector &secondary_vars,
            DecisionBuilder *db, const std::vector<SearchMonitor *> &monitors);

    /**
     * Prepare the solver before the call to Solve().
     *
     * Prepare involves setting the solver-validator connection, determining
     * local load (i.e., slices).
     *
     * @param validator validator to connect to
     */
    void Prepare(Validator &validator);

    /**
     * Creates an adapter to access the search array.
     *
     * @param name the adapater's name
     * @return access adapter
     */
    AdapterPtr CreateAdapter(const std::string &name);

    /**
     * Creates a new defaul SL search heuristic. The caller is responsible
     * for destruction.
     *
     * @param primary_vars primary decision variables
     * @param secondary_vars secondary variables
     *
     * @return SL search heuristic
     */
    DecisionBuilder *CreateDefaultHeuristic(const IntVarVector &primary_vars,
            const IntVarVector &secondary_vars);

    /**
     * Create default balancer for general heuristics.
     *
     * @param vars search variables
     * @param low low balance threshold
     * @param high high balance threshold
     * @return general balancer monitor
     */
    SearchMonitor *CreateBalancingMonitor(const IntVarVector &vars,
            double low, double high);

    /**
     * Returns user monitors attached to the main solver during the search.
     * Since this function returns a reference to the internal list, it is safe
     * to call before the Solve(). The list will be populated there. This also
     * means that checking the list before the search starts is meaningless,
     * since the list will be empty or incomplete.
     *
     * @return solver monitors defined by the user
     */
    const std::vector<SearchMonitor *> &GetUserMonitors() const {
        return search_monitors_.user_monitors_;
    }

    /**
     * Returns monitors established by Searchlight during the search.
     * Since this function returns a reference to the internal list, it is safe
     * to call before the Solve(). The list will be populated there. This also
     * means that checking the list before the search starts is meaningless,
     * since the list will be empty or incomplete.
     *
     * @return solver monitors defined by Searchlight
     */
    const std::vector<SearchMonitor *> &GetAuxMonitors() const {
        return search_monitors_.aux_monitors_;
    }

    /**
     * Returns the monitor that submits candidates to the validator. While it
     * is safe to call before Solve(), the function will return nullptr if
     * no monitor was established yet.
     *
     * @return the validator monitor or nullptr, if there is none
     */
    SearchMonitor *GetValidatorMonitor() const {
        return search_monitors_.validator_monitor_;
    }

    /**
     * Checks if any instances are available for help.
     *
     * @return true, if help available; false, otherwise
     */
    bool HelpAvailable() const {
        return !helpers_.empty();
    }

    /**
     * Handles a new helper for this solver.
     *
     * This handler also might reject help is the balancing is disabled (hard)
     * or if the solver is not prepared or searching (soft). The latter might
     * happen because of asynchronous messages, even when the solver has
     * already became idle.
     *
     * @param id helper  id
     */
    void HandleHelper(uint64_t id);

    /**
     * Dispatches work to an available helper.
     *
     * It is assumed that a helper is available, which can be checked by
     * the HelpAvailable() function.
     *
     * @param work assignments to off-load
     */
    void DispatchWork(const LiteAssignmentVector &work);

    /**
     * Prepares this solver as a helper.
     *
     * This function loads remote assignments into the solver.
     *
     * @param load remote load
     */
    void PrepareHelper(LiteAssignmentVector &load);

    /**
     * Return primary variables for the search.
     *
     * Primary variables are the variables that are used to divide the
     * search between instances, perform load balancing and make decisions
     * during the search.
     *
     * @return primary variables
     */
    const IntVarVector &GetPrimaryVars() const {
        return primary_vars_;
    }

    /**
     * Computes variable names. Should be used with care, since it parses
     * the entire model.
     */
    StringVector ComputeVarNames() const;

    /**
     * Connects validator to the solver.
     * @param validator vaidator to connect
     */
    void ConnectValidator(Validator *validator);

    /**
     * Main solver loop
     */
    void operator()();

    /**
     * Handles end-of-search, when all sovers are finished.
     */
    void HandleEndOfSearch() {
        std::lock_guard<std::mutex> lock{mtx_};
        status_ = Status::FINISHED;
        wait_cond_.notify_one();
        // We could join the thread here, but let's just postpone it to dtor
    }

    /**
     * Terminates the solver.
     */
    void HandleTerminate() {
        std::lock_guard<std::mutex> lock{mtx_};
        status_ = Status::TERMINATED;
        wait_cond_.notify_one();
    }

    /**
     * Starts the solver loop.
     */
    void Start() {
        thr_ = std::thread{std::ref(*this)};
    }

    /**
     * Return search (CP) solver.
     *
     * @return search solver
     */
    const Solver &GetSearchSolver() const {
        return solver_;
    }

    /**
     * Return search (CP) solver.
     *
     * @return search solver
     */
    Solver &GetSearchSolver() {
        return solver_;
    }

    /**
     * Return solver's status.
     *
     * @return solver's status
     */
    Status GetSolverStatus() const {
        return status_;
    }

private:
    // Reject/release helpers
    void RejectHelpers(bool hard);

    // Determines the local workload for this solver
    void DetermineLocalWorkload();

    // Check if need to terminate
    bool CheckTerminate() const {
        return status_ == Status::TERMINATED;
    }

    // Monitors participating in the search
    struct SearchMonitors {
        // Validator monitor for the search (owned by the solver -- no delete)
        ValidatorMonitor *validator_monitor_ = nullptr;

        // User monitors
        std::vector<SearchMonitor *> user_monitors_;

        // Auxiliary monitors established by SL
        std::vector<SearchMonitor *> aux_monitors_;

        // Returns a vector of all monitors that we install for search
        std::vector<SearchMonitor *> GetSearchMonitors() const {
            std::vector<SearchMonitor *> mons{user_monitors_};
            mons.insert(mons.end(), aux_monitors_.begin(), aux_monitors_.end());
            if (validator_monitor_) {
                mons.push_back(validator_monitor_);
            }
            return mons;
        }
    };

    // Searchlight
    Searchlight &sl_;

    // The solver
    Solver solver_;

    // Search heuristic
    DecisionBuilder *db_;

    // "Primary" decision variables -- defining the object
    IntVarVector primary_vars_;

    // "Secondary" decision variables -- conveying additional info
    IntVarVector secondary_vars_;

    // All variables, for convenience
    IntVarVector all_vars_;

    // Contains all vars: for prototyping and sending to validator
    Assignment vars_leaf_;

    // Monitors defined for the search
    SearchMonitors search_monitors_;

    // Total time spent on solving
    std::chrono::microseconds total_solve_time_{0};

    // Searchlight status
    std::atomic<Status> status_;

    // Helpers currently dispatched for this solver
    std::list<uint64_t> helpers_;

    // Work given to us by another instance
    LiteAssignmentVector helper_load_;

    // Local work defined at the beginning
    LiteAssignmentVector local_load_;

    // Adapters issued by the solver
    std::vector<AdapterPtr> adapters_;

    // If true, no local/remote assignment will be restored in Solve()
    bool local_remote_override_ = false;

    // Do we accept help?
    bool solver_balancing_enabled_ = true;

    // Solver id (top 32b: instance, low 32b: ordinal)
    uint64_t id_;

    // For loop control
    std::condition_variable wait_cond_;

    // Concurrency control
    std::mutex mtx_;

    // Solvers thread
    std::thread thr_;
};

/**
 * This class allows the search process to access data both via sampling
 * and real data. This class also provides the tools necessary to make
 * this access as efficient as possible. It provides a number of register
 * API functions via which the user can register search primitives. The rest
 * is handled by the Searchlight itself.
 */
class Searchlight : private boost::noncopyable {
public:
    /**
     * Status of Searchlight
     */
    enum class Status {
        ACTIVE,    //!< ACTIVE Solvers are active
        FIN_SEARCH,//!< FIN_SEARCH Solvers finished completely
        FIN_VALID, //!< FIN_VALID Validator finished local work
        COMMITTED, //!< COMMITTED SL completed successfully
        TERMINATED //!< TERMINATED Abnormal termination on error
    };

    /**
     * Maps UDF names to UDF creators.
     */
    typedef std::map<std::string, UDFFunctionCreator> UDFMapper;

    /**
     * Creates the main searchlight class. An instance of this class
     * corresponds to a single search process.
     *
     * @param task Searchlight task performing the query
     * @param dll_handler handler for loading/unloading DLLs
     */
    Searchlight(SearchlightTask &task,
            DLLHandler &dll_handler) :
        validator_(nullptr),
        validator_thread_(nullptr),
        array_desc_(nullptr),
        dl_udf_handle_(dll_handler.LoadDLL("searchlight_udfs")),
        sl_task_(task),
        status_(Status::ACTIVE) {}

    /**
     * Registers a data array and the corresponding sample for the search.
     *
     * We do not check the correspondence of the array and the sample, since
     * there are probably no means to do that.
     *
     * @param data the data array
     * @param samples samples available for the data array
     */
    void RegisterArray(ArrayPtr &data, const ArrayPtrVector &samples) {
        array_desc_.reset(new SearchArrayDesc(data, samples, *this));
    }

    /**
     * Registers an attribute for the search. All further adapter data
     * accesses must go through the returned id.
     *
     * This function also initializes all synopses corresponding to the
     * attribute.
     *
     * @param name the attribute's name
     * @return the access id for the attribute
     */
    AttributeID RegisterAttribute(const std::string &name);

    /**
     * Read and register query sequence from the specified file.
     *
     * @param search attribute id the sequence belongs to
     * @param filename file name to parse
     * @return registered query sequence id
     */
    size_t RegisterQuerySequence(AttributeID sattr, const std::string &filename);

    /**
     * Add tracking expression to Searchlight.
     *
     * Tracking expression's values are added to all solutions and output to
     * the users.
     *
     * The expression will be cast to var and the var will be assigned the
     * specified name. The overhead of adding an additional cast var to the
     * model should be negligent.
     *
     * @param expr tracking expression
     * @param name name to output
     */
    void AddTrackExpr(IntExpr *expr, const std::string &name);

    /**
     * Return vector of track expression names.
     *
     * The ordering in the vector might change if the new tracking vars are
     * added.
     *
     * @return track expression names
     */
    StringVector GetTrackExprs() const {
    	StringVector res;
    	std::copy(track_var_names_.begin(), track_var_names_.end(),
    			std::back_inserter(res));
    	return res;
    }

    /**
     * Return previously registered query sequence.
     *
     * The user must have previously registered the sequence and obtained
     * an id for it. This id is specified as the parameter.
     *
     * @param seq_id sequence id
     * @return query sequence reference
     */
    const DoubleVector &GetQuerySequence(size_t seq_id) const {
    	return query_seqs_.at(seq_id);
    }

    /**
     * Returns the creator for the requested UDF.
     *
     * @param name the NAME of the UDF
     * @return pointer to the creator function
     */
    UDFFunctionCreator GetUDFFunctionCreator(const std::string &name) {
        // first, look in the map
        std::string tag_name = "UDF_" + name;
        UDFFunctionCreator udf = GetRegisteredUDF(tag_name);
        if (udf) {
            return udf;
        }

        // else, look in the library
        std::string func_name = "Create" + tag_name;
        udf = (UDFFunctionCreator)dlsym(dl_udf_handle_,
                func_name.c_str());
        // We should check via dlerror, but NULL checking is fine
        if (!udf) {
            std::ostringstream err_msg;
            err_msg << "Cannot find a SL UDF function, name=" << func_name;
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << err_msg.str();
        }
        udf_map_[tag_name] = udf;
        return udf;
    }

    /**
     * Returns the creator for the registered UDF.
     *
     * @param tag_name the tag name for the UDF
     * @return the creator if found; NULL otherwise
     */
    UDFFunctionCreator GetRegisteredUDF(const std::string &tag_name) const {
        UDFMapper::const_iterator udf_it = udf_map_.find(tag_name);
        if (udf_it != udf_map_.end()) {
            return udf_it->second;
        } else {
            return NULL;
        }
    }

    /**
     * Return all UDFs requested from this SL.
     *
     * @return a set of UDF names requested from this SL.
     */
    StringSet GetAllUsedUDFs() const {
        StringSet res;
        for (UDFMapper::const_iterator cit = udf_map_.begin();
                cit != udf_map_.end(); cit++) {
            res.insert(cit->first);
        }
        return res;
    }

    /**
     * Causes the searchlight search to terminate. Note, it does not terminate
     * immediately, but within a reasonable amount of time.
     */
    void Terminate() {
        status_ = Status::TERMINATED;
        for (auto &s: solvers_) {
            s->HandleTerminate();
        }
    }

    /**
     * Checks if the searchlight is terminating.
     *
     * @return true, if we are terminating; false otherwise
     */
    bool CheckTerminate() const {
        return status_ == Status::TERMINATED;
    }

    /**
     * Returns the validator for the current search, if any.
     *
     * @return current validator or nullptr, if no exists
     */
    Validator& GetValidator() {
        return *validator_;
    }

    /**
     * Returns current Searchlight status.
     *
     * @return current Searchlight status
     */
    Status GetStatus() const {
        return status_;
    }

    /**
     * Reports locally finished validator.
     */
    void ReportFinValidator();

    /**
     * Handles the end-of-search request from the coordinator.
     *
     * This request basically means that all sorver processes finished. So,
     * we just change status and wake up validator if needed.
     */
    void HandleEndOfSearch();

    /**
     * Handles the commit request from the coordinator.
     *
     * This request basically means that everything is finished. We have to
     * change the status. Strictly speaking, we don't need to wake up the
     * validator -- it will be done at dtor. But we do it to reclaim
     * resources early.
     */
    void HandleCommit();

    /**
     * Converts the specified solution to a string representation.
     *
     * @param vals variable values
     * @return string representation of the solution
     */
    std::string SolutionToString(const std::vector<int64_t> &vals) const;

    /**
     * Returns config for the current task.
     *
     * @return task config
     */
    const SearchlightConfig &GetConfig() const;

    /**
     * Prepares Searchlight for execution.
     *
     * This method initiates every solver by using the task function,
     * starts the validator.
     *
     * @param name model name
     * @param task_fun user-defined task function
     * @param start_id starting solver id
     * @param solvers number of solvers to prepare
     */
    void Prepare(const std::string &name, SLTaskFunc task_fun,
            uint64_t start_id, int solvers);

    /**
     * Starts the solvers.
     */
    void StartSolvers() {
        for (auto &s: solvers_) {
            s->Start();
        }
    }

    /**
     * Handle a new helper for the solver.
     *
     * @param id helper id
     * @param solver local solver id (ordinal)
     */
    void HandleHelper(uint64_t id, uint32_t solver) {
        assert(solver < solvers_.size());
        solvers_[solver]->HandleHelper(id);
    }

    /**
     * Return specified solver's status.
     *
     * @param solver solver to get the status for
     * @return solver's status
     */
    SearchlightSolver::Status GetSolverStatus(uint32_t solver) const {
        assert(solver < solvers_.size());
        return solvers_[solver]->GetSolverStatus();
    }

    /**
     * Loads the helper-solver.
     *
     * @param load remote load
     * @param solver helper id
     */
    void PrepareHelper(LiteAssignmentVector &load, uint32_t solver) {
        assert(solver < solvers_.size());
        if (status_ != Status::TERMINATED) {
            solvers_[solver]->PrepareHelper(load);
        }
    }

    /**
     * Creates an adapter to access the search array.
     *
     * @param name the adapater's name
     * @return access adapter
     */
    AdapterPtr CreateAdapter(const std::string &name) {
        return std::make_shared<Adapter>(*array_desc_, name);
    }

    /**
     * Return search (CP) solver for the specified SL solver.
     *
     * @return search solver
     */
    const Solver &GetSearchSolver(uint32_t id) const {
        assert(id < solvers_.size());
        return solvers_[id]->GetSearchSolver();
    }

    /**
     * Return Searchligh solver with the given local id.
     *
     * primarily used to establish the model.
     *
     * @param id SL solver id
     * @return SL solver with the specified id
     */
    SearchlightSolver &GetSLSolver(uint32_t id) {
        assert(id < solvers_.size());
        return *solvers_[id];
    }

    /**
     * Cleanup Searchlight.
     *
     * This function should be called before the destructor.
     */
    void Cleanup();

private:

    friend class SearchlightSolver;

    // Solvers on this instance
    std::vector<std::unique_ptr<SearchlightSolver>> solvers_;

    // Var names
    StringVector var_names_;

    // Sequences registered with Searchlight (e.g., query waveform sequence)
    std::unordered_map<std::string, size_t> query_seqs_map_;
    std::vector<DoubleVector> query_seqs_;

    // Tracking variables
    std::unordered_set<std::string> track_var_names_;

    // Validator for the search and its thread
    Validator *validator_;
    std::thread *validator_thread_;

    // The array descriptor
    std::unique_ptr<SearchArrayDesc> array_desc_;

    // The udfs library
    void *dl_udf_handle_;

    // Maps requested UDF names to corresponding creators
    UDFMapper udf_map_;

    // Searchlight task
    SearchlightTask &sl_task_;

    // Searchlight status
    std::atomic<Status> status_;

    // For concurrency control
    std::mutex mtx_;
};

/**
 * Creates a cumulative time limit, applied to all nested searched combined.
 *
 * Users of the SL heuristic should use this time limit, since SL heuristic
 * might run a lot of nested searches. Unfortunately, or-tools doesn't have
 * an interface function that allows to set the ccumulativity. Thus, this
 * function.
 *
 * Users don't need to take care of freeing the memory. The solver will do it.
 *
 * @param s the solver to use the timer for
 * @param time_ms the cumulative time limit in ms
 * @return the time limit monitor
 */
SearchMonitor *MakeCumulativeTimeLimit(Solver &s, int64 time_ms);

} /* namespace searchlight */
#endif /* SEARCHLIGHT_SEARCHLIGHT_H_ */

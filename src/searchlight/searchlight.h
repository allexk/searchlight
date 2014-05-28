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

namespace searchlight {

class Validator;
class SearchlightCollector;

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
class DLLHandler {
public:
    /**
     * Constructs a new DLL handler. It assumes the default DLL directory
     * to be the SciDb plugin's directory.
     */
    DLLHandler() :
        dlls_dir_(scidb::Config::getInstance()->
                getOption<std::string>(scidb::CONFIG_PLUGINS)) {}

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
            std::string full_name = dlls_dir_ + "/lib" + name + ".so";
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
 * The type for a UDF function creator. It produces an or-tools IntExpr
 * representing the function. Takes as parameters: the solver to use
 * with, the adapter for accessing data, a vector of variables to work
 * with and integer parameters.
 */
typedef IntExpr *(* UDFFunctionCreator)(Solver *, AdapterPtr,
        const std::vector<IntVar *> &, const std::vector<int64> &);

/**
 * This class allows the search process to access data both via sampling
 * and real data. This class also provides the tools necessary to make
 * this access as efficient as possible. It provides a number of register
 * API functions via which the user can register search primitives. The rest
 * is handled by the Searchlight itself.
 */
class Searchlight {
public:
    /**
     * Maps UDF names to UDF creators.
     */
    typedef std::map<std::string, UDFFunctionCreator> UDFMapper;

    /**
     * Creates the main searchlight class. An instance of this class
     * corresponds to a single search process.
     *
     * @param name the name of the search
     */
    explicit Searchlight(const std::string &name, DLLHandler &dll_handler) :
        solver_(name),
        collector_(NULL),
        array_desc_(NULL),
        dl_udf_handle_(dll_handler.LoadDLL("searchlight_udfs")),
        terminate_(false) {}

    /**
     * The destructor.
     */
    ~Searchlight() {
        delete array_desc_;
    }

    /**
     * Returns the solver to use for the search.
     *
     * @return the solver for the search
     */
    Solver &GetSolver() {
        return solver_;
    }

    /**
     * Returns a constant reference to the solver.
     *
     * @return constant solver reference
     */
    const Solver &GetSolver() const {
        return solver_;
    }

    /**
     * Registers a data array and the corresponding sample for the search.
     *
     * We do not check the correspondence of the array and the sample, since
     * there are probably no means to do that.
     *
     * @param data the data array
     * @param sample the sample for the data array
     */
    void RegisterArray(ArrayPtr &data, ArrayPtr &sample) {
        array_desc_ = new SearchArrayDesc(data, sample);
    }

    /**
     * Registers an attribute for the search. All further adapter data
     * accesses must go through the returned id.
     *
     * @param name the attribute's name
     * @return the access id for the attribute
     */
    AttributeID RegisterAttribute(const std::string &name) {
        if (!array_desc_) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "No array registered with SL to register an attribute";
        }
        return array_desc_->RegisterAttribute(name);
    }

    /**
     * The main solve method that starts the search.
     *
     * @param db the decision builder (search heuristic)
     * @param vars the decision variables
     * @param monitors monitors, if required (can be empty)
     * @return true, if the search found something; false otherwise
     */
    bool Solve(DecisionBuilder *db, const IntVarVector &vars,
            const std::vector<SearchMonitor *> &monitors);

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
     * Creates an adapter to access the search array.
     *
     * @return access adapter
     */
    AdapterPtr CreateAdapter() const {
        return boost::make_shared<Adapter>(*array_desc_);
    }

    /**
     * Registers a solution collector for handling the exact results. SL is
     * not responsible for deleting it.
     *
     * @param collector collector for exact results
     */
    void RegisterCollector(SearchlightCollector *collector) {
        collector_ = collector;
    }

    /**
     * Causes the searchlight search to terminate. Note, it does not terminate
     * immediately, but within a reasonable amount of time.
     */
    void Terminate() {
        terminate_ = true;
    }

    /**
     * Checks if the searchlight is terminating.
     *
     * @return true, if we are terminating; false otherwise
     */
    bool CheckTerminate() const {
        /*
         * This code is obviously not thread-safe (see Terminate() as well),
         * but for our purposes it is fine. One thread will set it to true and
         * the other will read it. Delays in propagating the value are not that
         * important here, and it is changed only once.
         */
        return terminate_;
    }

private:
    // The solver
    Solver solver_;

    // Solution collector for main (exact) results
    SearchlightCollector *collector_;

    // The array descriptor
    SearchArrayDesc *array_desc_;

    // The udfs library
    void *dl_udf_handle_;

    // Maps requested UDF names to corresponding creators
    UDFMapper udf_map_;

    // True if we are required to terminate
    volatile bool terminate_;
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
     * @param validator the validator for checking assignments
     * @param vars a vector of decision variables (externally managed)
     * @param solver the main solver
     */
    ValidatorMonitor(Validator &validator, const IntVarVector &vars,
            Solver *solver) :
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

private:
    // The validator to pass the solution to
    Validator &validator_;

    // The vector of vars (managed outside)
    const IntVarVector &vars_;
};
} /* namespace searchlight */
#endif /* SEARCHLIGHT_SEARCHLIGHT_H_ */

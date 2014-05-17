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
#include "system/Config.h"

#include <dlfcn.h>
#include <boost/thread.hpp>
#include <boost/make_shared.hpp>

namespace searchlight {

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
    Searchlight(const std::string &name) :
        solver_(name),
        collector_(NULL),
        array_desc_(NULL),
        validator_(NULL),
        validator_thread_(NULL) {

        // loading the udf library
        const std::string &plugins_dir = scidb::Config::getInstance()->
                getOption<std::string>(scidb::CONFIG_PLUGINS);
        std::string lib_name = plugins_dir + "/searchlight_udfs.so";
        dl_udf_handle_ = dlopen(lib_name.c_str(), RTLD_LAZY | RTLD_LOCAL);
        if (!dl_udf_handle_) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Cannot load the UDF Searchlight library!";
        }
    }

    /**
     * The destructor.
     */
    ~Searchlight() {
        delete array_desc_;
        dlclose(dl_udf_handle_); // cannot be NULL

        delete validator_thread_;
        delete validator_;
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
    void RegisterArray(const Array &data, const Array &sample) {
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
    bool Solve(DecisionBuilder &db, const IntVarVector &vars,
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
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Cannot find a SL UDF function, name=" << func_name;
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
    void RegisterCollector(SolutionCollector *collector) {
        collector_ = collector;
    }

private:
    // The solver
    Solver solver_;

    // Solution collector for main (exact) results
    SolutionCollector *collector_;

    // The array descriptor
    SearchArrayDesc *array_desc_;

    // The udfs library
    void *dl_udf_handle_;

    // Maps requested UDF names to corresponding creators
    UDFMapper udf_map_;

    // The assignment validator and its thread
    Validator *validator_;
    boost::thread *validator_thread_;
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

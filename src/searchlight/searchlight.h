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
 * This is the main entry point for the search process. The Searchlight class
 * is responsible for registering all arrays, attributes, solvers and helping
 * the search process by providing the tools necessary.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_SEARCHLIGHT_H_
#define SEARCHLIGHT_SEARCHLIGHT_H_

#include "ortools_inc.h"
#include "array_desc.h"

namespace searchlight {

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
     * Creates the main searchlight class. An instance of this class
     * corresponds to a single search process.
     *
     * @param solver the or-tools main solver process
     */
    Searchlight(Solver &solver) :
        solver_(solver),
        array_desc_(NULL),
        array_adapter_(NULL) {}

    /**
     * The destructor.
     */
    ~Searchlight() {
        delete array_adapter_;
        delete array_desc_;
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
        array_adapter_ = new Adapter(*array_desc_);
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
     * Returns the search array's adapter.
     *
     * @return the adapter for the search array
     */
    const Adapter &GetArrayAdapter() const {
        return *array_adapter_;
    }

private:
    // The solver
    Solver &solver_;

    // The array descriptor
    SearchArrayDesc *array_desc_;

    // The array adapter
    const Adapter *array_adapter_;
};
} /* namespace searchlight */
#endif /* SEARCHLIGHT_SEARCHLIGHT_H_ */

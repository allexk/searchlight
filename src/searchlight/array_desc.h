/*
 * Copyright 2013, Brown University, Providence, RI.
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
 * @file array_desc.h
 * This file contains an array descriptor. It stores information about
 * dimensions, attributes, samples and real data for the given array.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_ARRAY_DESC_H_
#define SEARCHLIGHT_ARRAY_DESC_H_

#include "base.h"
#include "scidb_inc.h"
#include "array_access.h"

namespace searchlight {

/**
 * This is a descriptor for the given array. It contains information about
 * dimensions (names and intervals), attributes and other schema info. It also
 * contains pointers to the samples for individual arrays and the real data
 * fetched from the DBMS.
 *
 * It also contains a number of useful utility functions to deal with SciDb
 * arrays.
 */
class SearchArrayDesc {
public:
    /**
     * Creates a search array descriptor.
     *
     * @param array the scidb's original array
     * @param sample the sample array
     */
    SearchArrayDesc(const Array &array, const Array &sample) :
        array_(array),
        sampler_(sample, array.getArrayDesc()),
        data_accessor_(array) {}

    /**
     * Registers a search attribute with the descriptor. The sample for this
     * attribute is also loaded. This function returns the access id, and
     * all future attribute accesses via other SL classes should use this id.
     *
     * If the attribute has already been registered, this function returns
     * the same access id as previously.
     *
     * @param attr_name the name of the attribute to register
     * @return the access id for the attribute
     */
    AttributeID RegisterAttribute(const std::string &attr_name);

    /**
     * Returns the id of the specified attribute in the given vector. Note,
     * the attribute's id and the vector's id are not necessarily the same.
     * The id returned corresponds to the array's descriptor.
     *
     * @param attrs a vector of attributes
     * @param name the name of the attribute
     * @param res the id of the attribute (out)
     * @return true, if the id is found; false otherwise
     */
    static bool FindAttributeId(const Attributes &attrs,
            const std::string &name, AttributeID &res) {
        for (size_t i = 0; i < attrs.size; i++) {
            if (attrs[i].getName() == name) {
                res = attrs[i].getId();
                return true;
            }
        }
        return false;
    }

private:
    // The data array
    const Array &array_;

    /*
     *  Maps attribute names to search attribute IDs (might be different
     *  from the original attribute IDs).
     */
    std::vector<AttributeID> search_orig_ids_;

    // Maps  attribute names to internal access ids
    std::map<std::string, AttributeID> attr_to_id_;

    // The sampler
    Sampler sampler_;

    // The data accessor
    ArrayAccess &data_accessor_;
};

} /* namespace searchlight */
#endif /* SEARCHLIGHT_ARRAY_DESC_H_ */

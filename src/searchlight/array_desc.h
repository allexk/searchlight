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

#include "searchlight/base.h"
#include "searchlight/scidb_inc.h"

namespace searchlight {

/**
 * This is a descriptor for the given array. It contains information about
 * dimensions (names and intervals), attributes and other schema info. It also
 * contains pointers to the samples for individual arrays and the real data
 * fetched from the DBMS.
 */
class SearchArrayDesc {
public:
    /**
     * Returns the id of the specified attribute in the given vector. Note,
     * the attribute's id and the vector's id are not the same. The id
     * returned corresponds to the metadata.
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
};

} /* namespace searchlight */
#endif /* SEARCHLIGHT_ARRAY_DESC_H_ */

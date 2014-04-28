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
 * @file array_access.h
 * A class that allows accessing a SciDb data array.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_ARRAY_ACCESS_H_
#define SEARCHLIGHT_ARRAY_ACCESS_H_

#include "scidb_inc.h"

namespace searchlight {

/**
 * This class allows users to access SciDb arrays. For example, fetching
 * an element by the coordinates or an entire interval/region. It also allows
 * to iterate via elements of a region or call a callback on them. An
 * instance of this class corresponds to a single SciDb array.
 */
class ArrayAccess {
public:
    /**
     * Creates an accessor for a SciDb array.
     *
     * @param array the SciDb data array
     */
    ArrayAccess(const Array &array) : data_array_(array) {}
private:
    // The data array
    const Array &data_array_;
};
} /* namespace searchlight */
#endif /* SEARCHLIGHT_ARRAY_ACCESS_H_ */

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
 * @file adapter.h
 * A class that represent the main access point to the search array.
 * Depending on the situation (i.e., the mode of work) it might be
 * beneficial to access either the sample or the data itself. This class
 * does exactly that.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_ADAPTER_H_
#define SEARCHLIGHT_ADAPTER_H_

#include "scidb_inc.h"
#include "array_desc.h"

namespace searchlight {

/**
 * This class allows users to access search data. For example, fetching
 * an element by the coordinates or an entire interval/region. It also allows
 * to iterate via elements of a region or call a callback on them. Note, that
 * "data" might either correspond to a sample or the real array. The type
 * of access is determined by the adapter's mode. This class is exposed to
 * the library users for all search data access purposes.
 *
 * Note that data access functions return intervals, instead of single
 * values. This is done to make user's lives easier. An exact value is just
 * an interval of length 0. For example, asking for a(i, j) might give you
 * an exact value from the array, an approximate value generated from a sample
 * or an exact interval that contains a(i, j). This depends on the mode
 * of work, which can be selected explicitly, but is chosen automatically
 * by default depending on the situation.
 */
class Adapter {
public:
    /**
     * This is the mode of work of the adapter, which determines what
     * results it returns.
     */
    enum class Mode {
        /** The adapter gives the exact results from the real data */
        EXACT,
        /** The adapter gives approximate values generated from a sample */
        APPROX,
        /** The adapter gives the exact interval for the values */
        INTERVAL,
    };

    /**
     * Creates an adapter for a SciDb array.
     *
     * @param array the SciDb data array
     */
    Adapter(const SearchArrayDesc &array) :
        array_desc_(array),
        mode_(Mode::INTERVAL) {}

    /**
     * Sets the mode of operation for the adapter.
     *
     * @param mode the mode required
     */
    void SetAdapterMode(Mode mode) {
        mode_ = mode;
    }
private:
    // The search array descriptor
    const SearchArrayDesc &array_desc_;

    // Mode of operation
    Mode mode_;
};
} /* namespace searchlight */
#endif /* SEARCHLIGHT_ADAPTER_H_ */

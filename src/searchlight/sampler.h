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
 * @file smapler.h
 * This is a sampler. It allows access to an array sample stored in memory.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_SAMPLER_H_
#define SEARCHLIGHT_SAMPLER_H_

#include "searchlight/base.h"
#include "searchlight/array_desc.h"

#include "searchlight/scidb_inc.h"

#include <vector>

namespace searchlight {

/**
 * Sampler allows access to an array sample stored in memory. A sample is
 * a set of regions (chunks, strata), corresponding to different parts of the
 * array. The user access the sample using usual coordinates; the mapping
 * is achieved automatically. In general, sample regions might have an
 * arbitrary configuration, so the array intervals are no necessarily
 * evenly divided. To store the configuration in such a case an R-tree is
 * employed.
 *
 * Sampler corresponds only to a single attribute of the array. If more than
 * one attributes is uses, several Sampler instances is required.
 *
 * TODO: Go from the even case to the R-tree case.
 * TODO: We assume even regions covering the whole array without holes.
 */
class Sampler {
public:
    /**
     * Creates a sampler and load sample chunks from an array.
     * The array is  supposed to have two dimensions:
     *   first -- [0, number_of_chunks),
     *   second -- [0, original_attribute_ids). It also
     * The array is supposed to have two attributes:
     *   first -- min value for the chunk
     *   second -- max value for the chunk
     *
     * @param array the sample array
     */
    Sampler(const Array &array);

private:
    /*
     *  Parses chunk sizes out of the string. The string is suppposed to
     *  have the format "x_size,y_size".
     */
    void SetChunkSizes(const std::string &size_param);

    // A region of the sample
    struct Chunk {
        int64_t min_; // The minimum value in the region
        int64_t max_; // The maximum value in the region

        Chunk(int64_t min, int64_t max) : min_(min), max_(max) {}
    };

    // Sample chunks as a linearized array (multiple attributes)
    std::vector<std::vector<Chunk> > sample_chunks_;

    // The size of a chunk (one per dimension)
    int chunk_sizes_[2];

    Coordinates sample_start_;
};
} /* namespace searchlight */
#endif /* SEARCHLIGHT_SAMPLER_H_ */

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
 * @file sampler.cpp
 * The implementation of the sampler.
 *
 * @author Alexander Kalinin
 */

#include "sampler.h"

#include "array_desc.h"

namespace searchlight {

Sampler::Sampler(const Array &array) {
    const ArrayDesc &sample_desc = array.getArrayDesc();

    // by convenience we store sizes in the comment :)
    const std::string &sample_config = sample_desc.getComment();
    SetChunkSizes(sample_config);
    if (chunk_sizes_[0] == -1 || chunk_sizes_[1] == -1) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Cannot determine sample chunk sizes: sample="
                << sample_desc.getName();
    }

    // now, retrieve sample chunks for all attributes
    const Attributes &attrs = sample_desc.getAttributes(true);
    AttributeID min_id, max_id;
    if (!SearchArrayDesc::FindAttributeId(attrs, std::string("min"), min_id) ||
         !SearchArrayDesc::FindAttributeId(attrs, std::string("max"), max_id)) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Cannot find min/max attribute in the sample: sample="
                << sample_desc.getName();
    }
    boost::shared_ptr<ConstItemIterator> min_iterator =
            array.getItemIterator(min_id);
    boost::shared_ptr<ConstItemIterator> max_iterator =
            array.getItemIterator(max_id);

    // Sample: first dimension -- region, second -- the original attribute
    Coordinate chunk_end = sample_desc.getDimensions()[0].getCurrEnd();
    Coordinate attr_end = sample_desc.getDimensions()[1].getCurrEnd();
    Coordinates pos(2);
    sample_chunks_.resize(chunk_end);
    for (pos[0] = 0; pos[0] < chunk_end; pos[0]++) {
        for (pos[1] = 0; pos[1] < attr_end; pos[1]++) {
            min_iterator->setPosition(pos);
            max_iterator->setPosition(pos);
            const int64_t minv = min_iterator->getItem().getInt64();
            const int64_t maxv = max_iterator->getItem().getInt64();
            sample_chunks_[pos[0]].push_back(Chunk(minv, maxv));
        }
    }
}

void Sampler::SetChunkSizes(const std::string &size_param) {
    chunk_sizes_[0] = chunk_sizes_[1] = -1;
    sscanf(size_param.c_str(), "%d,%d", &chunk_sizes_[0], &chunk_sizes_[1]);
}

} /* namespace searchlight */

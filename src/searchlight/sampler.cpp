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

#include <boost/lexical_cast.hpp>

namespace searchlight {

Sampler::Sampler(const Array &array, const ArrayDesc &data_desc) :
        sample_array_(array),
        chunk_sizes_(data_desc.getDimensions().size()),
        sample_start_(data_desc.getDimensions().size()){
    const ArrayDesc &sample_desc = array.getArrayDesc();

    // by convenience we store sizes in the comment :)
    const std::string &sample_config = sample_desc.getComment();
    SetChunkSizes(sample_config);

    // The start of the sample corresponds to the start of the data array
    for (size_t i = 0; i < data_desc.getDimensions().size(); i++) {
        sample_start_[i] = data_desc.getDimensions()[i].getCurrStart();
    }

    // find min/max ids (not necessary, unless the array has the empty bitmap)
    const Attributes &attrs = sample_desc.getAttributes(true);
    if (!SearchArrayDesc::FindAttributeId(attrs, std::string("min"), min_id_) ||
        !SearchArrayDesc::FindAttributeId(attrs, std::string("max"), max_id_)) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Cannot find min/max attribute in the sample: sample="
                << sample_desc.getName();
    }

    if (sample_desc.getDimensions()[0].getCurrStart() != 0) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Chunk coordinate should start from 0: sample="
                << sample_desc.getName();
    }
    chunks_num_ = sample_desc.getDimensions()[0].getCurrEnd();
}

void Sampler::LoadSampleForAttribute(AttributeID attr_orig_id,
        AttributeID attr_search_id) {

    boost::shared_ptr<ConstItemIterator> min_iterator =
            sample_array_.getItemIterator(min_id_);
    boost::shared_ptr<ConstItemIterator> max_iterator =
            sample_array_.getItemIterator(max_id_);

    // Sample: first dimension -- region, second -- the original attribute
    sample_chunks_.push_back(ChunkVector());
    if (sample_chunks_.size() != attr_search_id) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Sampler and descriptor inconsistency: sample aid="
                << sample_chunks_.size() << ", desc id=" << attr_search_id;
    }
    ChunkVector &chunks = sample_chunks_.back();
    chunks.reserve(chunks_num_);
    Coordinates pos(2);
    pos[1] = attr_orig_id;
    for (pos[0] = 0; pos[0] < chunks_num_; pos[0]++) {
        min_iterator->setPosition(pos);
        max_iterator->setPosition(pos);
        const int64_t minv = min_iterator->getItem().getInt64();
        const int64_t maxv = max_iterator->getItem().getInt64();
        chunks.push_back(Chunk(minv, maxv));
    }
}

void Sampler::SetChunkSizes(const std::string &size_param) {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
    boost::char_separator<char> sep(",| ");
    tokenizer_t tokenizer(size_param, sep);

    int i = 0;
    for (tokenizer_t::const_iterator cit = tokenizer.begin();
            cit != tokenizer.end(); cit++) {
        chunk_sizes_[i++] = boost::lexical_cast<Coordinate>(cit->c_str());
    }
}
} /* namespace searchlight */

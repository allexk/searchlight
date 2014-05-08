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
 * @file array_access.cpp
 * The implementation of machinery for accessing SciDb arrays.
 *
 * @author Alexander Kalinin
 */

#include "array_access.h"

#include <query/ops/between/BetweenArray.h>
using scidb::BetweenArray;

using scidb::CountingAggregate;

namespace searchlight {

ValueVector ArrayAccess::ComputeAggreagate(const Coordinates &low,
        const Coordinates &high, AttributeID attr,
        const StringVector &aggr_names) const {
    // first, resolve the aggregates -- we use SciDb built-in ones
    SmallAggrVector aggrs;
    ValueVector res(aggr_names.size());
    bool need_nulls = false;
    for (size_t i = 0; i < aggr_names.size(); i++) {
        AggregatePtr agg = AggregateLibrary::getInstance()->
                createAggregate(aggr_names[i],
                        TypeLibrary::getType(attrs_[attr].getType()));
        aggrs.push_back(SmallAggr(agg, agg->isCounting(), !agg->ignoreNulls(),
                res[i]));
        if (!agg->ignoreNulls()) {
            need_nulls = true;
        }
    }

    // then, we use the SciDb's between array for the low-high region
    BetweenArray region(addEmptyTagAttribute(array_desc_), low, high,
            data_array_, tile_mode_ /* tile mode */);

    if (tile_mode_) {
        ComputeGeneralAggregateTile(region, attr, aggrs, need_nulls);
    } else {
        ComputeGeneralAggregate(region, attr, aggrs, need_nulls);
    }

    return res;
}

void ArrayAccess::ComputeGeneralAggregateTile(const Array &array,
        AttributeID attr, const SmallAggrVector &aggrs, bool need_nulls) {
    boost::shared_ptr<ConstArrayIterator> array_iter =
            array.getConstIterator(attr);

    int chunk_iter_mode = ConstChunkIterator::TILE_MODE |
            ConstChunkIterator::IGNORE_OVERLAPS |
            ConstChunkIterator::IGNORE_EMPTY_CELLS;
    if (!need_nulls) {
        chunk_iter_mode |= ConstChunkIterator::IGNORE_NULL_VALUES;
    }

    while (!array_iter->end()) {
        const ConstChunk &chunk = array_iter->getChunk();
        boost::shared_ptr<ConstChunkIterator> chunk_iter =
                chunk.getConstIterator(chunk_iter_mode);
        while (!chunk_iter->end()) {
            const Value &v = chunk_iter->getItem();
            const RLEPayload *tile = v.getTile();
            if (tile->count()) {
                for (size_t i = 0; i < aggrs.size(); i++) {
                    const SmallAggr &agg = aggrs[i];
                    if (agg.state_.getMissingReason() == 0) {
                        agg.agg_->initializeState(agg.state_);
                    }
                    agg.agg_->accumulatePayload(agg.state_, tile);
                }
            }
            ++(*chunk_iter);
        }
        ++(*array_iter);
    }
}

void ArrayAccess::ComputeGeneralAggregate(const Array &array,
        AttributeID attr, const SmallAggrVector &aggrs, bool need_nulls) {
    boost::shared_ptr<ConstArrayIterator> array_iter =
            array.getConstIterator(attr);

    int chunk_iter_mode = ConstChunkIterator::IGNORE_OVERLAPS |
            ConstChunkIterator::IGNORE_EMPTY_CELLS;
    if (!need_nulls) {
        chunk_iter_mode |= ConstChunkIterator::IGNORE_NULL_VALUES;
    }

    uint64_t total_count = 0, nonnull_count = 0;
    while (!array_iter->end()) {
        const ConstChunk &chunk = array_iter->getChunk();
        boost::shared_ptr<ConstChunkIterator> chunk_iter =
                chunk.getConstIterator(chunk_iter_mode);
        while (!chunk_iter->end()) {
            const Value &v = chunk_iter->getItem();
            total_count++;
            if (!v.isNull()) {
                nonnull_count++;
            }
            for (size_t i = 0; i < aggrs.size(); i++) {
                const SmallAggr &agg = aggrs[i];
                if (!agg.is_count_ && (agg.needs_nulls_ ||
                        !v.isNull())) {
                    if (agg.state_.getMissingReason() == 0) {
                        agg.agg_->initializeState(agg.state_);
                    }
                    agg.agg_->accumulate(agg.state_, v);
                }
            }
            ++(*chunk_iter);
        }
        ++(*array_iter);
    }

    // need to set count() if any
    for (size_t i = 0; i < aggrs.size(); i++) {
        const SmallAggr &agg = aggrs[i];
        if (agg.is_count_) {
            agg.agg_->initializeState(agg.state_);
            const uint64_t count = agg.needs_nulls_ ? total_count :
                    nonnull_count;
            ((CountingAggregate*)agg.agg_.get())->
                    overrideCount(agg.state_, count);
        }
    }
}

} /* namespace searchlight */

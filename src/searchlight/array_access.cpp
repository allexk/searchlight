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

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.array_access"));

ArrayAccess::~ArrayAccess() {
    const auto secs = std::chrono::duration_cast<std::chrono::duration<double>>(
            total_agg_time_).count();
    LOG4CXX_INFO(logger, "Total aggregates CPU time: " <<
            secs << 's');
}

TypedValueVector ArrayAccess::ComputeAggreagate(const Coordinates &low,
        const Coordinates &high, AttributeID attr,
        const StringVector &aggr_names) const {
    // first, resolve the aggregates -- we use SciDb built-in ones
    SmallAggrVector aggrs;
    TypedValueVector res(aggr_names.size());
    bool need_nulls = false;
    for (size_t i = 0; i < aggr_names.size(); i++) {
        AggregatePtr agg = AggregateLibrary::getInstance()->
                createAggregate(aggr_names[i],
                        TypeLibrary::getType(attrs_[attr].getType()));
        aggrs.push_back(SmallAggr(agg, agg->isCounting(), !agg->ignoreNulls()));
        if (!agg->ignoreNulls()) {
            need_nulls = true;
        }
    }

    /*
     * Then, we use the SciDb's between array for the low-high region.
     *
     * Caveat: BetweenArray uses ++ somewhat frivolously. If we have
     * a DepartArray as the input (surely for distributed environment),
     * it might cause a number of round-trips to other instances. Data
     * won't be fetched; only inquires about some remote chunks.
     *
     * We're not using SubArray here (also an appropriate choice), since it
     * doesn't support tile iteration. Performance wise, it should be around
     * the same, although it uses iterators more carefully.
     *
     * FIXME: Rewrite BetweenArray to use setPosition() instead of ++
     */
    scidb::SpatialRangesPtr between_ranges =
    		boost::make_shared<scidb::SpatialRanges>(low.size());
    between_ranges->_ranges.push_back(scidb::SpatialRange(low, high));
    const BetweenArray region(array_desc_, between_ranges, data_array_);

    if (tile_mode_) {
        ComputeGeneralAggregateTile(region, attr, aggrs, need_nulls);
    } else {
        ComputeGeneralAggregate(region, attr, aggrs, need_nulls);
    }

    // finalize aggregates
    for (size_t i = 0; i < aggrs.size(); i++) {
        const SmallAggr &agg = aggrs[i];
        agg.agg_->finalResult(res[i].first, agg.state_);
        res[i].second = agg.agg_->getResultType();
    }

    return res;
}

void ArrayAccess::ComputeGeneralAggregateTile(const Array &array,
        AttributeID attr, SmallAggrVector &aggrs, bool need_nulls) const {
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

        // starting the timer
        const auto agg_start_time = std::chrono::steady_clock::now();

        while (!chunk_iter->end()) {
            const Value &v = chunk_iter->getItem();
            const RLEPayload *tile = v.getTile();
            if (tile->count()) {
                for (size_t i = 0; i < aggrs.size(); i++) {
                    SmallAggr &agg = aggrs[i];
                    agg.agg_->accumulateIfNeeded(agg.state_, tile);
                }
            }
            ++(*chunk_iter);
        }

        // stopping the timer
        const auto agg_end_time = std::chrono::steady_clock::now();

        total_agg_time_ += std::chrono::duration_cast<decltype(total_agg_time_)>(
                agg_end_time - agg_start_time);

        ++(*array_iter);
    }
}

void ArrayAccess::ComputeGeneralAggregate(const Array &array,
        AttributeID attr, SmallAggrVector &aggrs, bool need_nulls) const {
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

        // starting the timer
        const auto agg_start_time = std::chrono::steady_clock::now();

        while (!chunk_iter->end()) {
            const Value &v = chunk_iter->getItem();
            total_count++;
            if (!v.isNull()) {
                nonnull_count++;
            }
            for (size_t i = 0; i < aggrs.size(); i++) {
                SmallAggr &agg = aggrs[i];
                if (!agg.is_count_ && (agg.needs_nulls_ ||
                        !v.isNull())) {
                    agg.agg_->accumulateIfNeeded(agg.state_, v);
                }
            }
            ++(*chunk_iter);
        }

        // stopping the timer
        const auto agg_end_time = std::chrono::steady_clock::now();

        total_agg_time_ += std::chrono::duration_cast<decltype(total_agg_time_)>(
                agg_end_time - agg_start_time);

        ++(*array_iter);
    }

    // need to set count() if any
    for (size_t i = 0; i < aggrs.size(); i++) {
        SmallAggr &agg = aggrs[i];
        if (agg.is_count_) {
            agg.agg_->initializeState(agg.state_);
            const uint64_t count = agg.needs_nulls_ ? total_count :
                    nonnull_count;
            ((CountingAggregate*)agg.agg_.get())->
                    overrideCount(agg.state_, count);
        }
    }
}

bool ArrayAccess::SqDistance(const Coordinates &low, const Coordinates &high,
		size_t int_dim, AttributeID attr,
		const DoubleVector &query_seq, double &res) const {
	// Interval coordinates
	const Coordinate low_int = low[int_dim];
	const Coordinate high_int = high[int_dim];
	if (low_int > high_int || (high_int - low_int + 1) != query_seq.size()) {
		return false;
	}
    // Iterate through the interval
	ConstItemIteratorPtr iter = data_array_->getItemIterator(attr, 0);
    const Type &value_type = TypeLibrary::getType(attrs_[attr].getType());
    res = 0;
    Coordinates pos{low};
    for (Coordinate i = low_int; i <= high_int; ++i) {
    	pos[int_dim] = i;
    	if (!iter->setPosition(pos)) {
    		// Empty element
    		return false;
    	}
		Value &value = iter->getItem();
		if (value.isNull()) {
			return false;
		}
		// Non-empty and non-NULL
		const double dval = scidb::ValueToDouble(value_type.typeId(), value);
    	const double val_diff = query_seq[i - low_int] - dval;
    	res += val_diff * val_diff;
    }
    return true;
}
} /* namespace searchlight */

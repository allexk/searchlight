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
    explicit ArrayAccess(ArrayPtr array) :
        data_array_(array),
        array_desc_(array->getArrayDesc()),
        attrs_(array_desc_.getAttributes(false)),
        last_attr_(-1) {

        /*
         * We use tiles only if the storage uses them, since we work with SciDb
         * arrays (borrowed from the optimizer).
         */
        tile_mode_ = Config::getInstance()->
                        getOption<bool>(scidb::CONFIG_RLE_CHUNK_FORMAT) &&
                     Config::getInstance()->
                        getOption<int>(scidb::CONFIG_TILE_SIZE) > 1;
    }

    /**
     * Computes aggregates over a region of the array. The aggregate is taken
     * from the SciDb aggregates library, based on its name and the
     * attribute's type. count(*) is also supported for including NULL elements.
     * In the latter case, any AttributeID can be specified.
     *
     * @param low the leftmost corner of the region
     * @param high the rightmost corner of the region
     * @param attr the id of the attribute (SciDb one)
     * @param aggr_names names of the aggregates
     * @return results, one per corresponding aggregate
     */
    TypedValueVector ComputeAggreagate(const Coordinates &low,
            const Coordinates &high, AttributeID attr,
            const StringVector &aggr_names) const;

    /**
     * Returns the element at the given position. The element itself is
     * returned in the corresponding parameter, and the function itself
     * returns the status of access. Out-of-bounds accesses are treated as
     * empty elements. A NULL value might be returned, and it is the
     * responsibility of the caller to handle that.
     *
     * @param point element's coordinates
     * @param attr requested attribute
     * @param res reference to assign the element to
     * @return true, if the element is returned;
     *         false, if empty or out-of-bounds
     */
    bool GetElement(const Coordinates &point, AttributeID attr,
            TypedValue &res) const {
        const ConstItemIteratorPtr &array_iter = last_iter_;
        if (attr != last_attr_) {
            ItemIteratorsMap::const_iterator it = item_iters_.find(attr);
            if (it != item_iters_.end()) {
                array_iter = it->second;
                last_attr_ = attr;
            } else {
                ConstItemIteratorPtr iter = data_array_->getItemIterator(attr,
                        ConstChunkIterator::IGNORE_OVERLAPS |
                        ConstChunkIterator::IGNORE_EMPTY_CELLS);
                array_iter = item_iters_[attr] = iter;
                last_attr_ = attr;
            }
        }

        if (array_iter->setPosition(point)) {
            // non-empty guaranteed, but can be NULL
            res.first = array_iter->getItem();
            res.second = TypeLibrary::getType(attrs_[attr].getType());
            return true;
        }
        return false; // empty or out-of-bounds
    }

private:
    /*
     *  The structure represents an internal "small" aggregate. This is done
     *  basically to handle multi-aggregate requests.
     */
    struct SmallAggr {
        AggregatePtr agg_;        // the aggregate itself
        bool is_count_;           // is it count()/count(*)
        bool needs_nulls_;        // does it need NULLs?
        Value state_;            // the accumulator -- NULL(0) by default

        SmallAggr(AggregatePtr agg, bool is_count, bool needs_nulls) :
            agg_(agg),
            is_count_(is_count),
            needs_nulls_(needs_nulls) {}
    };
    typedef std::vector<SmallAggr> SmallAggrVector;

    // Aggregate computation for tile-based inputs.
    static void ComputeGeneralAggregateTile(const Array &array,
            AttributeID attr, const SmallAggrVector &aggrs, bool need_nulls);

    // General aggregate for non-tile computations
    static void ComputeGeneralAggregate(const Array &array, AttributeID attr,
            const SmallAggrVector &aggrs, bool need_nulls);

    // The data array
    ArrayPtr data_array_;

    // The array's descriptor
    const ArrayDesc &array_desc_;

    // The array's attributes
    const Attributes &attrs_;

    // Do we use tiles for aggregates?
    bool tile_mode_;

    // Iterators for point access
    typedef std::map<AttributeID, ConstItemIteratorPtr> ItemIteratorsMap;
    ItemIteratorsMap item_iters_;

    // Lats used iterator for faster access
    ConstItemIteratorPtr last_iter_;
    AttributeID last_attr_;
};
} /* namespace searchlight */
#endif /* SEARCHLIGHT_ARRAY_ACCESS_H_ */

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
#include "sampler.h"

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
     * Structure describing chunk zones for the distribution.
     *
     * Each zone defines successive divisions of the original array into
     * slices along a single coordinate. Thus, two numbers: start (origin) and
     * the slice size.
     */
    struct ChunkZones {
        struct Zone {
            Coordinate start_;
            uint64_t slice_;
            size_t dim_num_;
        };
        std::vector<Zone> zones_;
    };

    /**
     * Creates a search array descriptor. We pass a vector of samples to it.
     * The first in the vector becomes the primary sample, while others
     * become auxiliary samples to optimize region reads.
     *
     * @param array the scidb's original array
     * @param samples available samples
     * @param sl_config Searchlight config
     */
    SearchArrayDesc(const ArrayPtr &array, const ArrayPtrVector &samples,
            const SearchlightConfig &sl_config) :
        array_(array),
        sampler_(array->getArrayDesc(), samples, sl_config),
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
     * @param load_aux_samples load auxiliary samples for the attribute
     * @return the access id for the attribute
     */
    AttributeID RegisterAttribute(const std::string &attr_name,
            bool load_aux_samples);

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
        for (size_t i = 0; i < attrs.size(); i++) {
            if (attrs[i].getName() == name) {
                res = attrs[i].getId();
                return true;
            }
        }
        return false;
    }

    /**
     * Returns a sampler for this array.
     *
     * @return this array's sampler
     */
    const Sampler &GetSampler() const {
        return sampler_;
    }

    /**
     * Returns the data accessor for this array.
     *
     * @return this array's data accessor
     */
    const ArrayAccess &GetAccessor() const {
        return data_accessor_;
    }

    /**
     * Returns original (i.e., SciDb) array descriptor.
     *
     * @return original array descriptor
     */
    const ArrayDesc &GetOriginalArrayDesc() const {
        return array_->getArrayDesc();
    }

    /**
     * Returns the original array attribute id for the specified internal
     * access id.
     *
     * @param attr the internal attribute id
     * @return the original attribute id
     */
    AttributeID GetArrayAttrributeID(AttributeID attr) const {
        return search_orig_ids_[attr];
    }

    /**
     * Fills in chunk distribution info for a set of chunks.
     *
     * Chunk distribution is information about which instance has which
     * chunks. It is passed to the function as a vector of counters (one per
     * instance), and the vector is updated by using information about
     * static distribution and chunks fetched by the local instance during
     * the query execution.
     *
     * @param query caller's query context
     * @param chunk_pos set of chunk positions
     * @param distr vector of instance counters (will be updated)
     */
    void GetDynamicChunksDistribution(const boost::shared_ptr<Query> &query,
            const CoordinateSet &chunk_pos, std::vector<int> &distr) const;

    /**
     * Fills in stripe-based chunk distribution for a set of chunks.
     *
     * Stripe distribution gives every validator a stripe of chunks, which
     * is determined by splitting the first array dimension. It is assumed
     * that the number of stripes is equal to the number of counters in the
     * input array.
     *
     * By default, this functions uses dimension from the zone specification
     * to select the appropriate coordinate for each of the chunk positions.
     * However, if coordinates are one-dimensional, the dimension is ignored
     * and the only available coordinate is used instead. This is especially
     * helpful in case of transferring forwards, when it is inefficient to
     * send full coordinates.
     *
     * The function can work with any sequential container that supports
     * forward traversing. It is assumed that the container contains
     * coordinates in vector form.
     *
     * @param chunk_pos set of chunk positions
     * @param distr output counters (one per active instance)
     * @param zone distribution to use
     */
    template <typename CoordinatesSequence>
    void GetStripesChunkDistribution(const CoordinatesSequence &chunk_pos,
            std::vector<int> &distr, const ChunkZones::Zone &zone) const;

    /**
     * Create chunk zones according to the specified parameters.
     *
     * The function starts from the whole array area. Then, at step i,
     * it picks the maximum length interval and slices the
     * dimension into slice_nums[i] pieces. This describes a zone.
     * If more steps is needed, it contracts the current area by
     * picking the slice with number slice_ords[i].
     *
     * The number of zones is determined by the size of slice_nums.
     *
     * @param slice_nums number of slices at each step
     * @param slice_ords the slice to pick for the next iteration
     * @return zones
     */
    ChunkZones CreateChunkZones(const std::vector<size_t> &slice_nums,
            const std::vector<size_t> &slice_ords) const;

    /**
     * Return the entire search array size.
     *
     * This function takes only registered attributes into account, not all
     * array attributes.
     *
     * @return search array size
     */
    size_t GetSearchArraySize() const;

private:
    // The data array
    const ArrayPtr array_;

    /*
     *  Maps internal access ids to original attribute IDs.
     */
    std::vector<AttributeID> search_orig_ids_;

    // Maps  attribute names to internal access ids
    std::unordered_map<std::string, AttributeID> attr_to_id_;

    // The sampler
    Sampler sampler_;

    // The data accessor
    ArrayAccess data_accessor_;
};

} /* namespace searchlight */
#endif /* SEARCHLIGHT_ARRAY_DESC_H_ */

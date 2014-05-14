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
 * @file sampler.h
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
     * Creates a sampler and loads the info about the sample. The info
     * includes the number of chunks, the starting point, chunk
     * attribute ids.
     *
     * The array is  supposed to have two dimensions:
     *   first -- [0, number_of_chunks),
     *   second -- [0, original_attribute_ids). It also
     * The array is supposed to have attributes:
     *   first  -- min value for the chunk
     *   second -- max value for the chunk
     *   third  -- density of the chunk
     *
     * @param array the sample array
     * @param data_desc the descriptor of the data array
     */
    Sampler(const Array &array, const ArrayDesc &data_desc);

    /**
     * Loads the sample for the particular attribute. It is assumed that the
     * sample array's attribute ids correspond to the data array ones.
     *
     * @param attr_orig_id the attribute id as of the data array
     * @param attr_search_id the internal attribute id for the search
     */
    void LoadSampleForAttribute(AttributeID attr_orig_id,
            AttributeID attr_search_id);


    /**
     * Registers a new sample aggregate.
     *
     * @param name aggregate's name
     * @param aggr aggregate's factory
     */
    void RegisterAggregate(const std::string &name,
            SampleAggregateFactory *aggr) {
        aggrs_[name] = aggr;
    }

    /**
     * Computes specified aggregates for the specified regions for the
     * specified attribute. The aggregates must be registered in the sampler or
     * be default ones. The boundaries returned for the aggregate are guaranteed
     * to be precise and reachable, i.e., they are not over-/under-estimated.
     * As for the approximate value, it is computed by using the uniformity
     * assumption, which might change in the future.
     *
     * @param low the low coordinates of the region
     * @param high the upper coordinates for the region
     * @param attr the attribute to compute for
     * @param aggr_names the names of the aggregates
     * @return results, one per aggregate in the form of intervals; see the
     *  definition of the interval for details
     */
    IntervalValueVector ComputeAggregate(const Coordinates &low,
            const Coordinates &high, AttributeID attr,
            const StringVector &aggr_names) const;

    /**
     * Returns an approximate value at the specified coordinates. The value is
     * computed using the uniformity assumption. The min/max boundaries are
     * precise.
     *
     * @param point the point's coordinates
     * @param attr the attribute to compute the value for
     * @return the approximate value with the boundaries
     */
    IntervalValue GetElement(const Coordinates &point, AttributeID attr) const;

    /**
     * A region (chunk) of the sample.
     */
    struct Chunk {
        /**
         * The minimum value in the region
         */
        double min_;

        /**
         * The maximum value in the region.
         */
        double max_;

        /**
         * The sum of the elements in the region.
         */
        double sum_;

        /**
         * The number of non-empty/non-null elements.
         */
        uint64_t count_;

        /**
         * Constructs a new chunk.
         *
         * @param min the minimum value in the chunk
         * @param max the maximum value in the chunk
         * @param sum the sum of all values in the chunk
         * @param count the number of non-empty/non-null elements
         */
        Chunk(double min, double max, double sum, uint64_t count) :
            min_(min), max_(max), sum_(sum), count_(count) {}

        /**
         * Checks if the chunk is empty. Empty means all its elements are
         * empty or null.
         *
         * @return true, if the chunk is empty; false otherwise
         */
        bool empty() const {
            return count_ == 0;
        }
    };

    /**
     * A vector of sample chunks.
     */
    typedef std::vector<Chunk> ChunkVector;

private:
    /*
     *  Iterator over the chunks of a region. We assume that chunks are
     *  laid in the row-major order for the purpose of returning the
     *  linear number of the current chunk the iterator points to.
     */
    class RegionIterator {
    public:
        /*
         * Creates an iterator over the region, specified by low and
         * high coordinates (both inclusive).
         *
         * We assume that the low-high coordinates comprise a valid region,
         * where high[i] >= low[i].
         */
        RegionIterator(const Sampler &sampler, const Coordinates &low,
                const Coordinates &high) :
                    pos_(low),
                    chunk_pos_(-1),
                    region_low_(low),
                    region_high_(high),
                    sampler_(sampler),
                    valid_(true) {
            chunk_pos_ = GetChunkPos();
        }

        // prefix ++
        RegionIterator &operator++() {
            size_t i = pos_.size() - 1;
            while ((pos_[i] += sampler_.chunk_sizes_[i]) > region_high_[i]) {
                pos_[i] = region_low_[i];
                if (i == 0) {
                    valid_ = false;
                    break;
                }
                i--;
            }

            if (valid_) {
                if (i == pos_.size() - 1) {
                    chunk_pos_++;
                } else {
                    chunk_pos_ = GetChunkPos();
                }
            }

            return *this;
        }

        // Is the iterator invalid (at the end?)
        bool end() const {
            return !valid_;
        }

        // Return the linear chunk position
        uint64_t GetChunk() const {
            return chunk_pos_;
        }

        // Are we fully covering the current chunk with the region?
        bool CoversFullChunk() const {
            bool res = true;
            for (size_t i = 0; i < pos_.size(); i++) {
                if (pos_[i] % sampler_.chunk_sizes_[i] != 0) {
                    res = false;
                    break;
                }
            }
            return res;
        }

        /*
         *  Return the covered portion of the current chunk. If the
         *  chunk is fully covered it will return its coordinates.
         */
        void GetCoveredChunk(Coordinates &low, Coordinates &high) const {
            low = pos_;
            for (size_t i = 0; i < pos_.size(); i++) {
                const Coordinate chunk_size = sampler_.chunk_sizes_[i];
                high[i] = (pos_[i] - pos_[i] % chunk_size) + chunk_size - 1;
            }
        }

        /*
         * Computes the number of elements in the covered part of the chunk.
         */
        uint64_t GetPartSize() const {
            uint64_t res = 1;
            for (size_t i = 0; i < pos_.size(); i++) {
                const Coordinate chunk_size = sampler_.chunk_sizes_[i];
                res *= chunk_size - pos_[i] % chunk_size;
            }
            return res;
        }

    private:
        // Returns the linear chunk number for the current position
        uint64_t GetChunkPos() const {
            Coordinates chunk_coord(pos_);
            for (size_t i = 0; i < chunk_coord.size(); i++) {
                chunk_coord[i] = (chunk_coord[i] - sampler_.sample_origin_[i]) /
                        sampler_.chunk_sizes_[i];
            }

            // TODO: single out one-/two-dimensional cases?
            uint64_t pos = 0;
            for (size_t i = 0; i < chunk_coord.size(); i++) {
                pos *= sampler_.chunk_nums_[i];
                pos += chunk_coord[i];
            }

            return pos;
        }

        // Current position
        Coordinates pos_;

        // Current chunk position
        uint64_t chunk_pos_;

        // Boundaries for the region
        const Coordinates &region_low_, &region_high_;

        // The sampler we are traversing
        const Sampler &sampler_;

        // Do we point to a valid position?
        bool valid_;
    };

    /*
     *  Parses chunk sizes out of the string. The string is suppposed to
     *  have the format "x_size,y_size,...".
     */
    void SetChunkSizes(const std::string &size_param);

    /*
     * Checks that specified region bounds are correct (low <= high) and
     * aligns them appropriately, to fit within the sampled area.
     */
    void CheckAndCorrectBounds(Coordinates &low, Coordinates &high) const;

    /*
     * Checks if the point's coordinates are within the sample area.
     */
    bool CheckBounds(const Coordinates &point) const;

    /*
     * Computes the number of elements in a chunk based on its length.
     */
    uint64_t ComputeShapeChunkSize() const;

    // Attribute IDs for min/max/count/sum elements in the sample array
    AttributeID min_id_, max_id_, count_id_, sum_id_;;

    // The number of sample chunks
    Coordinate chunks_num_;

    // The sampple array
    const Array &sample_array_;

    // Sample chunks as a linearized array (multiple attributes)
    std::vector<ChunkVector> sample_chunks_;

    // The size of a chunk (one per dimension)
    Coordinates chunk_sizes_;

    // The number of elements in a chunk, based on its shape
    uint64_t shape_chunk_size_;

    // The starting and ending points of the sample
    Coordinates sample_origin_, sample_end_;

    // The number of chunks per each dimension
    Coordinates chunk_nums_;

    // Aggregate map to reso;ve aggregates
    typedef std::map<const std::string, SampleAggregateFactory> AggrMap;
    AggrMap aggrs_;
};

/**
 * This is a base class for all sample aggregates. Sample means it returns
 * an interval result in the form of [min, max] values and a possible
 * indication of NULL.
 */
class SampleAggregate {
public:
    /**
     * Adds another chunk to the aggregate. Since the aggregate is computed
     * for a region, that region might intersect only a part of the
     * corresponding sample chunk. To account for this, the function is
     * called with the number of elements in the intersection. Only non-empty
     * chunks calls this function. So, there might be a case when only the
     * finalize is called -- when the whole region is definitely empty.
     *
     * @param chunk_size the number of elements in the corresponding chunk
     * @param part_size the number of elements in the intersection
     * @param chunk reference to the corresponding sample chunk
     */
    virtual void AccumulateChunk(uint64_t chunk_size,
            uint64_t part_size, const Sampler::Chunk &chunk) = 0;

    /**
     * Finalizes the result after all the chunks have been traversed. For
     * example, this is the place to compute avg by dividing sum by count.
     *
     * @param res reference to put the final result
     */
    virtual void Finalize(IntervalValue &res) = 0;

private:
    /**
     * Destructor.
     */
    virtual ~SampleAggregate() {}
};

/**
 * Sample Aggregate Factory.
 */
typedef SampleAggregate *(*SampleAggregateFactory)();

/**
 * A vectort of shared pointers to sample aggregates.
 */
typedef std::vector<boost::shared_ptr<SampleAggregate>> SampleAggregatePtrVector;


} /* namespace searchlight */
#endif /* SEARCHLIGHT_SAMPLER_H_ */

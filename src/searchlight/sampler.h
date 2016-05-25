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
 * This is a sampler. It allows access to array synopses.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_SAMPLER_H_
#define SEARCHLIGHT_SAMPLER_H_

#include "base.h"
#include "scidb_inc.h"
#include "cache.h"
#include "synopsis.h"

#include <mutex>

namespace searchlight {

class SampleAggregate;

/**
 * A vector of shared pointers to sample aggregates.
 */
typedef std::vector<std::shared_ptr<SampleAggregate>> SampleAggregatePtrVector;

/**
 * Sample Aggregate Factory.
 */
typedef SampleAggregate *(*SampleAggregateFactory)();

/**
 * This class contains all synopsis information about a single attribute
 * of the array. This includes synopsis cells and the corresponding
 * metadata. Different synopses, even for a single attribute,
 * are represented by different Synopsis objects.
 */
class Synopsis : private boost::noncopyable,
                 public GridSynopsis<AggCell, AggCellItemReader> {
private:
    typedef GridSynopsis<AggCell, AggCellItemReader> Base;

public:
    /**
     * Create a new synopsis.
     *
     * No cells are loaded at this point. Only meta-data is initialized.
     *
     * @param data_desc data array descriptor
     * @param array synopsis array
     */
    Synopsis(const ArrayDesc &data_desc, const ArrayPtr &array);

    /**
     * Check if the region can be computed by this synopsis.
     *
     * @param low low region coordinates
     * @param high high region coordinates
     * @return true, if this synopses can compute the region;
     *   false, otherwise
     */
    bool RegionValidForSample(const Coordinates &low,
            const Coordinates &high) const;

    /**
     * Compute the region's aggregates over this synopsis.
     *
     * This function does not "finalyze" the aggregates. This is the
     * responsibility of the caller.
     *
     * @param low low region coordinates
     * @param high high region coordinates
     * @param aggs resolver aggregate functions
     */
    void ComputeAggregate(const Coordinates &low,
            const Coordinates &high,
            const SampleAggregatePtrVector &aggs);

    /**
     * Get the interval for array element.
     *
     * @param point element coordinates
     * @return element interval
     */
    IntervalValue GetElement(const Coordinates &point);

    /**
     * Compute an aggregate for regions considering cell threshold.
     *
     * If a region intersects a cell in less than the threshold part
     * (specified as a [0, 1] ratio), the intersection is put into the
     * specified list and is not computed. All other intersections
     * (above the threshold) are computed as usual.
     *
     * Note, if the threshold is <=0, this synopsis will compute all
     * regions entirely, without leftovers.
     *
     * @param aggs aggregates to compute
     * @param in_regions input regions to compute
     * @param left_regions regions left without the computation
     * @param cell_thr cell threshold
     */
    void ComputeAggregatesWithThr(const SampleAggregatePtrVector &aggs,
            const std::vector<Region> &in_regions,
            std::vector<Region> &left_regions, double cell_thr);
};
using SynopsisPtr = std::unique_ptr<Synopsis>;
using SynopsisSharedPtr = std::shared_ptr<Synopsis>;

/**
 * This class contains DFT-based synopsis information about 1D
 * sequences. Its cells are actually MBRs of the
 * corresponding DFT trace subsequences.
 *
 * This type of synopsis has a number of conventions:
 *
 * 1) Its name is <name>_<dft_params>_<usual_chunk_size>
 *      Where <dft_params_is>: <subseq_sizeXdft_componentsXMBR_size>
 * 2) The last dimension stores DFT coordinates (i.e., its not grid)
 * 3) The dimension before last is the waveform one
 */
class DFTSynopsis : private boost::noncopyable,
                    public GridSynopsis<DFTCell, DFTCellItemReader> {
private:
    typedef GridSynopsis<DFTCell, DFTCellItemReader> Base;

public:
    /**
     * Create a new DFT synopsis.
     *
     * No cells are loaded at this point. Only meta-data is initialized.
     *
     * @param data_desc data array descriptor
     * @param array synopsis array
     */
    DFTSynopsis(const ArrayDesc &data_desc, const ArrayPtr &array);

    /**
     * Returns the MBR size of the synopsis.
     *
     * The MBR size is equivalent to the trace interval used to construct
     * the synopsis MBRs.
     *
     * @return synopsis MBR size
     */
    uint64_t GetMBRSize() const {
        return mbr_size_;
    }

    /**
     * Compute minimum square distance from the point to interval MBRs.
     *
     * Look at all MBRs covering the interval and compute minimal
     * square distance to the query point.
     *
     * @param low low interval boundary
     * @param high high interval boundary
     * @param point query points
     * @return minumum square distance from the interval to the point
     */
    IntervalValue SqDist(Coordinate low, Coordinate high,
            const DoubleVector &points);

    /**
     * Return subsequence size (omega) for this synopsis.
     * @return synopsis subsequence size
     */
    size_t GetSubsequenceSize() const {
        return subseq_size_;
    }

    /**
     * Return the number of DFT coordinates in the synopsis.
     *
     * @return number of DFT coordinates
     */
    size_t GetDFTNum() const {
        return dft_num_;
    }

private:
    /*
     * Memory required to store each cell.
     *
     * Overriden for DFT, since it has additional storage.
     */
    virtual size_t MemoryPerCell() const override {
        // In general, we store 2 vectors for MBRs
        return sizeof(DFTCell) + sizeof(double) * dft_num_ * 2;
    }

    // Check if the point inside the managed interval
    void CheckBounds(Coordinate &point) const;

    // Parse DFT-synopsis specific params
    void ParseDFTSynopsisParams(const std::string &params);

    // Subsequence size for which DFT coordinates were computed
    size_t subseq_size_;

    // The number of DFT coordinates.
    uint64_t dft_num_;

    // The size of the MBR (i.e., the trace size)
    size_t mbr_size_;
};
using DFTSynopsisPtr = std::unique_ptr<DFTSynopsis>;

/**
 * Sampler allows access to array synopses stored in memory. A synopsis is
 * a set of regions (cells), corresponding to different parts of the
 * array. The user access the synopsis using usual coordinates; the mapping
 * is achieved automatically. In general, synopsis regions might have an
 * arbitrary configuration, so the array intervals are no necessarily
 * evenly divided.
 *
 * TODO: Go from the even case to the R-tree case.
 * TODO: We assume even regions covering the whole array without holes.
 */
class Sampler {
private:
    /**
     * Type of caching for the synopsis.
     */
    enum class CachingType {
        EAGER, /** All cells are in a linear array, high footprint. */
        LAZY,  /** All _requested_ chunks are in a set of arrays. */
        NONE   /** No caching. Each cell is read via array iterators.*/
    };

public:
    /**
     * Creates a sampler and a catalog of all synopses available.
     *
     * The catalog is indexed by the attribute name. Each array's name is
     * supposed to consist of three components: <name>_<attr>_NxN... where
     * name is arbitrary, attr is the attribute name and Ns define synopsis
     * cell size.
     *
     * An array itself is supposed to be of the same dimensionality as the
     * original array, with indexes scaled down appropriately depending on
     * the cell size. Attributes of the array are: min, max, sum, count. The
     * order is not important.
     *
     * @param data_desc data array descriptor
     * @param synopsis_arrays synopsis arrays
     * @param sl_config Searchlight config
     */
    Sampler(const ArrayDesc &data_desc,
            const std::vector<ArrayPtr> &synopsis_arrays,
            const SearchlightConfig &sl_config);

    /**
     * Destructor.
     *
     * For now it just outputs some statistics.
     */
    ~Sampler();

    /**
     * Loads synopses for the particular attribute.
     *
     * If the synopses has already been loaded, the function does nothing.
     * If, on the other hand, the attribute is new, synopses are created
     * (meta-data only) and sorted in descending order of their cell sizes.
     * The most "coarse" is used as a "primary" synopsis for all
     * estimations.
     *
     * This function enables caching only for the primary synopsis. All others
     * will not cache any cells, unless explicitly told otherwise.
     *
     * @param attr_name the attribute name
     * @param attr_search_id the internal attribute id for the search
     */
    void LoadSampleForAttribute(const std::string &attr_name,
            AttributeID attr_search_id);

    /**
     * Clears the persistent synopsis cache.
     *
     * Note, all current clients using synopses previously cached won't be
     * affected. However, when they release their arrays, those will be lost
     * until created again. This is expected functionality.
     */
    static void ClearPersistentCache();

    /**
     * Registers a new sample aggregate.
     *
     * @param name aggregate's name
     * @param aggr aggregate's factory
     */
    void RegisterAggregate(const std::string &name,
            SampleAggregateFactory aggr) {
        aggrs_[name] = aggr;
    }

    /**
     * Register query sequence for the specified attribute.
     *
     * Basically we compute DFTs based on the synopses available and cache
     * the for future use.
     *
     * @param attr sequence attribute
     * @param seq_id sequence id
     * @param seq query sequence
     */
    void RegisterQuerySequence(AttributeID attr, size_t seq_id,
    		const DoubleVector &seq);

    /**
     * Computes specified aggregates for the specified region for the
     * specified attribute. The aggregates must be registered in the sampler or
     * be default ones. The boundaries returned for the aggregate are guaranteed
     * to be precise and reachable, i.e., they are not over-/under-estimated.
     * As for the approximate value, it is computed by using the uniformity
     * assumption, which might change in the future.
     *
     * This function can also try to compute the exact value of the
     *
     * @param low the low coordinates of the region
     * @param high the upper coordinates for the region
     * @param s_attr the search (internal) attribute
     * @param aggr_names the names of the aggregates
     * @return results, one per aggregate in the form of intervals; see the
     *  definition of the interval for details
     */
    IntervalValueVector ComputeAggregate(const Coordinates &low,
            const Coordinates &high, AttributeID s_attr,
            const StringVector &aggr_names, bool exact) const;

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
     * Compute minimum square distance from the point to interval MBRs.
     *
     * Look at all MBRs covering the interval and compute minimal
     * square distance to the query point.
     *
     * This function assumes the query sequence has been registered with the
     * sampler. The id must be provided as a parameter by the user.
     *
     * @param attr attribute to compute the distance for
     * @param low low interval boundary
     * @param high high interval boundary
     * @param seq_id query sequence id
     * @return minimum square distance from the interval to the point
     */
    IntervalValue SqDist(AttributeID attr, Coordinate low, Coordinate high,
    		size_t seq_id) const;

    /**
     * Return available DFT synopsis sizes.
     *
     * When the query sequence is transformed into a DFT point, it is
     * necessary to know which synopses the sampler has, since each synopsis is
     * created for the fixed (parameter) size. This function returns available
     * sizes.
     *
     * The sizes are returned only for the specified attribute id
     *
     * @param internal attribute id
     * @return vector of available DFT synopsis sizes
     */
    std::vector<size_t> AvailableDFTSizes(AttributeID attr) const {
    	std::vector<size_t> res;
    	if (attr < dft_synopses_.size() && !dft_synopses_[attr].empty()) {
        	// FIXME: For now we use only a single DFT synopsis
    		res.push_back(dft_synopses_[attr][0]->GetSubsequenceSize());
    	}
    	return res;
    }

private:
    // Synopsis cache
    using SynopsisCache = SearchlightCache<std::string, Synopsis>;

    /*
     * Every synopsis array name consists of:
     *   <user-defined name>_<attr_name>_NxN... where N are synopsis cell
     *   sizes.
     *
     *   This function just returns the three parts in separate strings.
     */
    static StringVector ParseArrayParamsFromName(const std::string &array_name);

    /*
     * Prepares synopses (sorts by cell size, caches, preloads).
     *
     * skip_synopses specifies synopses to SKIP from prepare. Others are
     * prepared, which is usefule for inter-query cached synopses.
     */
    template<class T>
    void PrepareSynopses(std::vector<T> &synopses,
            const std::unordered_set<size_t> &skip_synopses,
    		CachingType cache_type, bool preload, size_t mem_limit,
			const std::string &attr_name, AttributeID attr_search_id);

    /*
     * Compute DFTs of size dft_size from the sequence seq and store first
     * dft_num components in res sequentially.
     *
     * Since DFT results in complex numbers, we store real and imaginary
     * parts sequentially for components. For example, if dft_num == 6,
     * we take the first 3 components of the DFT.
     */
    static void ComputeDFTs(const DoubleVector &seq, size_t dft_size,
    		size_t dft_num, DoubleVector &res);

    // Synopses catalog (attr. name --> list of synopsis)
    std::unordered_map<std::string,std::vector<ArrayPtr>> array_synopses_;

    /*
     * Catalog of prepared synopses. First vector is indexed by internal
     * attribute ids, second vector contains synopses ranged by the cell
     * size in decreasing order.
     */
    std::vector<std::vector<SynopsisSharedPtr>> synopses_;
    std::vector<std::vector<DFTSynopsisPtr>> dft_synopses_;

    /*
     * Searchlight inter-query cache for synopses.
     */
    static SynopsisCache synopsis_cache_;

    // Descriptor of the data array we store synopses for
    const ArrayDesc data_desc_;

    // Searchlight config
    const SearchlightConfig &sl_config_;

    // Aggregate map to resolve aggregates
    typedef std::map<const std::string, SampleAggregateFactory> AggrMap;
    AggrMap aggrs_;

    // Cell threshold: if the query region covers a cell less, we go deeper
    double cell_thr_;

    // MBR threshold: if the query region is covered less, we'll go deeper
    double mbr_thr_;

    // Limit on the total number of cells to use for estimations
    size_t cell_limit_;

    // Do we cache synopses between queries?
    bool cache_synopses_;

    /*
     *  Query sequences DFT cache types and routines.
     */
    // DFTs are cached by attribute, sequence id and synopsis resolution
    using DFTSequenceID = std::tuple<AttributeID, size_t, size_t>;
    // Hash function to use in maps
    struct DFTSequenceIDHash : public std::unary_function<DFTSequenceID, size_t> {
    	size_t operator()(const DFTSequenceID &id) const {
    		size_t h = 17;
    		h = 31 * h + std::get<0>(id);
    		h = 31 * h + std::get<1>(id);
    		h = 31 * h + std::get<2>(id);
    		return h;
    	}
    };
    // Cache of DFT sequences
    using DFTSequenceCache = std::unordered_map<DFTSequenceID, DoubleVector, DFTSequenceIDHash>;
    // Cache itself
    DFTSequenceCache dft_seq_cache_;
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
            uint64_t part_size, const AggCell &chunk) = 0;

    /**
     * Finalizes the result after all the chunks have been traversed. For
     * example, this is the place to compute avg by dividing sum by count.
     *
     * @param res reference to put the final result
     */
    virtual void Finalize(IntervalValue &res) = 0;

    /**
     * Destructor.
     */
    virtual ~SampleAggregate() {}
};

} /* namespace searchlight */
#endif /* SEARCHLIGHT_SAMPLER_H_ */

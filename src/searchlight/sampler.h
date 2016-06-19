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
class AggSynopsis : private boost::noncopyable,
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
    AggSynopsis(const ArrayDesc &data_desc, const ArrayPtr &array);

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
using AggSynopsisSharedPtr = std::shared_ptr<AggSynopsis>;

/**
 * This class is a sequence-based synopsis (i.e., DFT/PAA type).
 * Its cells are actually MBRs of the corresponding trace subsequences.
 *
 * This type of synopsis has a number of conventions:
 *
 * 1) Its name is <name>_<params>_<usual_chunk_size>
 *      Where <params> is: <name(dft/paa)Xsubseq_size>
 * 2) The last dimension stores components (i.e., its not a grid one),
 *      from which we take the number of components
 * 3) The dimension before last is the waveform one. Thus <usual_chunk_size>
 *      contains the size of the MBR.
 */
class SeqSynopsis : private boost::noncopyable,
                    public GridSynopsis<SeqCell, SeqCellItemReader> {
private:
    typedef GridSynopsis<SeqCell, SeqCellItemReader> Base;

public:
    enum class Type {
        DFT,  // DFT synopsis
        PAA   // PAA synopsis
    };

    /**
     * Create a new sequence synopsis.
     *
     * No cells are loaded at this point. Only meta-data is initialized.
     *
     * @param data_desc data array descriptor
     * @param array synopsis array
     */
    SeqSynopsis(const ArrayDesc &data_desc, const ArrayPtr &array);

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
     * @param seq_info information about the query sequence
     * @return minumum square distance from the interval to the point
     */
    IntervalValue SqDist(const Coordinates &low, const Coordinates &high,
            const TransformedSequenceInfo &seq_info);

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
    size_t GetFeaturesNum() const {
        return features_num_;
    }

    /**
     * Return this synopsis type.
     *
     * @return this synopsis type
     */
    Type GetType() const {
        return type_;
    }

private:
    /*
     * Memory required to store each cell.
     *
     * Overriden for DFT, since it has additional storage.
     */
    virtual size_t MemoryPerCell() const override {
        // In general, we store 2 vectors for MBRs
        return sizeof(SeqCell) + sizeof(double) * features_num_ * 2;
    }

    // Subsequence size for which the features were computed
    size_t subseq_size_;

    // The number of features.
    uint64_t features_num_;

    // The size of the MBR (i.e., the trace size)
    size_t mbr_size_;

    // Sequence synopsis type
    Type type_;
};
using SeqSynopsisSharedPtr = std::shared_ptr<SeqSynopsis>;

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
     * @param low low interval boundary
     * @param high high interval boundary
     * @param seq_id query sequence id
     * @return minimum square distance from the interval to the point
     */
    IntervalValue SqDist(const Coordinates &low, const Coordinates &high,
    		size_t seq_id) const;

private:
    // Synopsis cache
    using AggSynopsisCache = SearchlightCache<std::string, AggSynopsis>;
    using SeqSynopsisCache = SearchlightCache<std::string, SeqSynopsis>;

    /*
     * Prepares synopses (sorts by cell size, caches, preloads).
     *
     * skip_synopses specifies synopses to SKIP from prepare. Others are
     * prepared, which is usefule for inter-query cached synopses.
     */
    template<class SynType>
    void PrepareSynopses(std::vector<SynType> &synopses,
    		CachingType cache_type, bool preload, size_t mem_limit,
			const std::string &attr_name, AttributeID attr_search_id) {
        // First, sort synopses by their cell sizes
        std::stable_sort(synopses.begin(), synopses.end(),
                [](const SynType &s1, const SynType& s2) {
            return s1->GetCellSize() > s2->GetCellSize();
        });

        // Then set cache type and preload
        LOG4CXX_INFO(logger_, "Caching synopses with memory limit " <<
                     mem_limit << "MB");
        for (size_t i = 0; i < synopses.size(); ++i) {
            const auto &syn = synopses[i];
            const size_t syn_mem_mb = syn->MemorySize() / 1024 / 1024; // in MB
            if (!mem_limit || syn_mem_mb <= mem_limit) {
                if (syn->GetCachingType() == CachingType::NONE) {
                    // Cache the synopsis
                    syn->SetCacheType(cache_type);
                    // .. and preload it if needed
                    if (preload) {
                        syn->Preload();
                    }
                } else {
                    LOG4CXX_WARN(logger_, "Synopsis " << syn->GetName() <<
                            " has already been cached. Skipping.");
                }
                if (mem_limit) {
                    mem_limit -= syn_mem_mb;
                }
            } else {
                LOG4CXX_WARN(logger_, "Cannot cache synopsis " <<
                        syn->GetName() << ". Memory is exhausted: left=" <<
                        mem_limit << "MB");
                break;
            }
        }

        // Warn about performance problems...
        if (!synopses.front()->IsCached()) {
            LOG4CXX_WARN(logger_, "No synopses are cached for attribute "
                    << attr_name << "(" << attr_search_id << ")");
        }

        // Debug printing
        if (logger_->isDebugEnabled()) {
            std::ostringstream msg;
            msg << "Synopses loaded for attribute " << attr_name
                    << '(' << attr_search_id << "): ";
            for (const auto &syn: synopses) {
                msg << syn->GetName() << "(cached: " << int(syn->GetCachingType())
                        << "), ";
            }
            logger_->debug(msg.str());
        }
    }

    /*
     * Compute DFTs of size dft_size from the sequence seq and store first
     * dft_num components in res sequentially.
     *
     * Since DFT results in complex numbers, we store real and imaginary
     * parts sequentially for components. For example, if dft_num == 6,
     * we take the first 3 components of the DFT.
     */
    static void ComputeDFTs(const DoubleVector &seq, size_t ss_size,
    		size_t dft_num, DoubleVector &res);

    /*
     * Compute PAA for the specified sequence.
     */
    static void ComputePAA(const DoubleVector &seq, size_t ss_size,
                           size_t feat_num, DoubleVector &res);

    // Synopses catalog (attr. name --> list of synopsis)
    std::unordered_map<std::string,std::vector<ArrayPtr>> array_synopses_;

    /*
     * Catalog of prepared synopses. First vector is indexed by internal
     * attribute ids, second vector contains synopses ranged by the cell
     * size in decreasing order.
     */
    struct AttributeSynopses {
        bool loaded_ = false;
        std::vector<AggSynopsisSharedPtr> agg_synopses_;
        std::vector<SeqSynopsisSharedPtr> seq_synopses_;
    };
    std::vector<AttributeSynopses> synopses_;

    /*
     * Searchlight inter-query cache for synopses.
     */
    static AggSynopsisCache agg_synopsis_cache_;
    static SeqSynopsisCache seq_synopsis_cache_;

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
     *  Query sequences cache types and routines.
     */
    using SequenceCache = std::unordered_map<size_t, TransformedSequenceInfo>;
    SequenceCache seq_cache_;

    // Logger
    log4cxx::LoggerPtr logger_;
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

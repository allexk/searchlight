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
 * @file adapter.h
 * A class that represent the main access point to the search array.
 * Depending on the situation (i.e., the mode of work) it might be
 * beneficial to access either the sample or the data itself. This class
 * does exactly that.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_ADAPTER_H_
#define SEARCHLIGHT_ADAPTER_H_

#include "scidb_inc.h"
#include "base.h"

#include <chrono>
#include <unordered_set>

namespace searchlight {

class SearchArrayDesc;

/**
 * Approximate value, returned by an estimator.
 *
 * TODO: there might be some benefits in returning the probability of MAY_NULL
 * in the future.
 */
struct IntervalValue {
    /**
     * Describes the state of the value
     */
    enum State {
        NON_NULL,//!< NON_NULL definitely not NULL
        MAY_NULL,//!< MIN_NULL might be NULL (thus, min is invalid)
        NUL      //!< NUL definitely NULL
    };

    /**
     * Minimum possible
     */
    double min_;

    /**
     * Maximum possible
     */
    double max_;

    /**
     * An approximate value.
     */
    double val_;

    /**
     * Value state
     */
    State state_;

    /**
     * Creates a NULL interval value.
     */
    IntervalValue() : min_(0), max_(0), val_(-1), state_(NUL) {}
};

/**
 * Outputs an IntervalValue into an output stream.
 *
 * @param os the stream to output to
 * @param iv the value to output
 * @return the stream used for the output
 */
std::ostream &operator<<(std::ostream &os, const IntervalValue &iv);

/**
 * A vector of interval values.
 */
typedef std::vector<IntervalValue> IntervalValueVector;

/**
 * This class allows users to access search data. For example, fetching
 * an element by the coordinates or an entire interval/region. It also allows
 * to iterate via elements of a region or call a callback on them. Note, that
 * "data" might either correspond to a sample or the real array. The type
 * of access is determined by the adapter's mode. This class is exposed to
 * the library users for all search data access purposes.
 *
 * Note that data access functions return intervals, instead of single
 * values. This is done to make user's lives easier. An exact value is just
 * an interval of length 0. For example, asking for a(i, j) might give you
 * an exact value from the array, an approximate value generated from a sample
 * or an exact interval that contains a(i, j). This depends on the mode
 * of work, which can be selected explicitly, but is chosen automatically
 * by default depending on the situation.
 */
class Adapter {
public:
    /**
     * This is the mode of work of the adapter, which determines what
     * results it returns.
     */
    enum Mode {
        /** The adapter gives the exact results from the real data */
        EXACT = 0,
        /** The adapter gives approximate values generated from a sample */
        APPROX,
        /** The adapter gives a precise interval of possible values */
        INTERVAL,
        /** The adapter always gives (-inf, +inf) estimations */
        DUMB,
    };

    /**
     * Statistics about data accesses. To make sense, collecting stats must
     * be enabled.
     */
    struct AccessStats {
        /**
         * Start positions of chunks accessed.
         */
        std::unordered_set<Coordinates, CoordinatesHash> chunks_pos_;

        /**
         * Reset statistics.
         */
        void Reset() {
            chunks_pos_.clear();
        }
    };

    /**
     * Creates an adapter for a SciDb array.
     *
     * @param array the SciDb data array
     */
    Adapter(const SearchArrayDesc &array, const std::string &name) :
        array_desc_(array),
        mode_(Mode::DUMB),
        name_(name),
        stats_enabled_(false) {}

    /**
     * Destructor.
     */
    ~Adapter();

    /**
     * Sets the mode of operation for the adapter.
     *
     * @param mode the mode required
     */
    void SetAdapterMode(Mode mode) {
        mode_ = mode;
    }

    /**
     * Computes specified aggregates for the specified regions for the
     * specified attribute. Depending on the mode, the result will be
     * approximate with precise boundaries or exact, computed over the real
     * data. For the exact result, the min/max/val of each result are the same.
     * All aggregates must be registered with the sampler and with SciDb.
     * Default ones are registered automatically.
     *
     * @param low the low coordinates of the region
     * @param high the upper coordinates for the region
     * @param attr the access id of the attribute to compute
     * @param aggr_names the names of the aggregates
     * @return results, one per aggregate in the form of intervals; see the
     *  definition of the interval for details
     */
    IntervalValueVector ComputeAggregate(const Coordinates &low,
            const Coordinates &high, AttributeID attr,
            const StringVector &aggr_names) const;

    /**
     * Returns the value at the specified coordinates. The value can be
     * approximate or exact based on this adapter's mode of operation. In
     * case of the exact value the min/max/val of the result is the same.
     *
     * @param point the point's coordinates
     * @param attr the attribute to compute the value for
     * @return the value
     */
    IntervalValue GetElement(const Coordinates &point, AttributeID attr) const;

    /**
     * Return square distance from the query sequence to the interval values.
     *
     * The value can be approximated by an interval via the DFT sampling with
     * MBRs, covering traces. The exact value is computed by retrieving array
     * values belonging to the interval.
     *
     * Note, that the query interval might be represented by several points,
     * which corresponds to dividing the query sequence into several
     * subsequence. The latter is needed to be able to use DFT sampling --
     * it can handle only sequences if the specified length.
     *
     * It is assumed that the query sequence has been registered with
     * Searchlight, so the sequence id is specified as a parameter.
     *
     * @param low left interval boundary
     * @param high right interval boundary
     * @param int_dim interval dimension among low/high
     * @param seq_id query sequence id
     * @return square distance from point to interval values
     */
    IntervalValue SqDist(const Coordinates &low, const Coordinates &high,
    		size_t int_dim, AttributeID attr, size_t seq_id) const;

    /**
     * Enables stats collecting mode.
     *
     * Statistics is collected about all accesses performed during the time
     * period the mode is active. The stats is also reset.
     *
     */
    void StartCollectingStats() {
        stats_enabled_ = true;
        stats_.Reset();
    }

    /**
     *
     * Disables collecting stats about accesses.
     *
     * The stats is not reset, so the caller might retrieve it via
     * subsequent calls.
     */
    void StopCollectingStats() {
        stats_enabled_ = false;
    }

    /**
     * Returns current statistics about data accesses.
     *
     * @return current stats
     */
    const AccessStats &GetCurrentStats() const {
        return stats_;
    }

    /**
     * Returns search array descriptor this adapter handles.
     *
     * @return search array descriptor of this adapter
     */
    const SearchArrayDesc &GetSearchArrayDesc() const {
        return array_desc_;
    }

private:
    /*
     * Iterates over chunks of the specified region. Chunks correspond to
     * the specified array. The iterator is initially set and then steps over
     * starting positions of the chunks.
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
        RegionIterator(const ArrayDesc &array_desc, const Coordinates &low,
                const Coordinates &high) :
                    array_desc_(array_desc),
                    pos_(low),
                    region_low_(low),
                    region_high_(high),
                    valid_(true) {
            // Correct coordinates
            const scidb::Dimensions &dims = array_desc_.getDimensions();
            for (size_t i = 0; i < dims.size(); i++) {
                const Coordinate low = dims[i].getCurrStart();
                const Coordinate high = dims[i].getCurrEnd();
                if (region_low_[i] < low) {
                    region_low_[i] = low;
                    pos_[i] = low;
                }
                if (region_high_[i] > high) {
                    region_high_[i] = high;
                }
            }
            // Align coords, since we're interested in start positions only
            array_desc_.getChunkPositionFor(pos_);
            array_desc_.getChunkPositionFor(region_low_);
            array_desc_.getChunkPositionFor(region_high_);
        }

        // prefix ++
        RegionIterator &operator++() {
            const scidb::Dimensions &dims = array_desc_.getDimensions();
            size_t i = pos_.size() - 1;
            while ((pos_[i] += dims[i].getChunkInterval()) > region_high_[i]) {
                pos_[i] = region_low_[i];
                if (i == 0) {
                    valid_ = false;
                    break;
                }
                i--;
            }
            return *this;
        }

        // Is the iterator invalid (at the end?)
        bool end() const {
            return !valid_;
        }

        // Returns current position
        const Coordinates &CurrentPosition() const {
            assert(valid_);
            return pos_;
        }

    private:
        // The original array descriptor
        const ArrayDesc &array_desc_;

        // Current position
        Coordinates pos_;

        // Boundaries for the region
        Coordinates region_low_, region_high_;

        // Are we pointing to a valid position?
        bool valid_;
    };

    void UpdateStatsWithRegion(const Coordinates &low,
            const Coordinates &high) const;

    // The search array descriptor
    const SearchArrayDesc &array_desc_;

    // Mode of operation
    Mode mode_;

    // The adapter's name
    const std::string name_;

    // Total time spent for requests
    mutable std::chrono::microseconds total_req_time_{0};

    // Accesses statistics collected in stats mode
    mutable AccessStats stats_;

    // True, if collecting stats is enabled
    bool stats_enabled_;
};

/**
 * Adapter shared pointer.
 */
typedef std::shared_ptr<Adapter> AdapterPtr;

} /* namespace searchlight */
#endif /* SEARCHLIGHT_ADAPTER_H_ */

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
#include "adapter.h"
#include "scidb_inc.h"

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
     * A region (cell) of the synopsis.
     */
    struct Cell {
        /**
         * The minimum value in the region
         */
        double min_ = 0;

        /**
         * The maximum value in the region.
         */
        double max_ = 0;

        /**
         * The sum of the elements in the region.
         */
        double sum_ = 0;

        /**
         * The number of non-empty/non-null elements.
         */
        uint64_t count_ = 0;

        /**
         * Cell has been loaded.
         */
        std::atomic<bool> valid_{false};

        /**
         * Constructs a new cell.
         *
         * @param min the minimum value in the chunk
         * @param max the maximum value in the chunk
         * @param sum the sum of all values in the chunk
         * @param count the number of non-empty/non-null elements
         */
        Cell(double min, double max, double sum, uint64_t count) :
            min_(min), max_(max), sum_(sum), count_(count), valid_{true} {}

        /**
         * Copy constructor.
         *
         * Needed because of the atomic and usage of this class in a vector.
         *
         * @param cell cell to copy from
         */
        Cell(const Cell &cell) :
            min_(cell.min_), max_(cell.max_), sum_(cell.sum_),
            count_(cell.count_),
            valid_(cell.valid_.load(std::memory_order_relaxed)) {}

        /**
         * Constructs an invalid (default) cell.
         *
         * If the cell is invalid, that means it has to be loaded from disk.
         */
        Cell() {}
    };

    /**
     * A vector of sample chunks.
     */
    typedef std::vector<Cell> Cells;

private:
    class RegionIterator;

    /**
     * This class contains all synopsis information about a single attribute
     * of the array. This includes synopsis cells and the corresponding
     * metadata. Different synopses, even for a single attribute,
     * are represented by different Synopsis objects.
     */
    class Synopsis : private boost::noncopyable {
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
         * Enable or disable caching of cells.
         *
         * @param mode true, enable caching; false, disable caching
         */
        void SetCacheMode(bool mode);

        /**
         * Returns the cell size of the synopsis.
         *
         * @return synopsis cell size
         */
        uint64_t GetCellSize() const {
            return shape_cell_size_;
        }

        /**
         * Check that specified region bounds are correct (low <= high) and
         * align them appropriately, to fit within the sampled area.
         *
         * @param low low coordinates of the region
         * @param high high coordinates of the region
         */
        void CheckAndCorrectBounds(Coordinates &low, Coordinates &high) const;

        /**
         * Check if the point's coordinates are within the synopsis area.
         *
         * @param point point's coordinates
         * @return true, if the point is within the area; false, otherwise
         */
        bool CheckBounds(const Coordinates &point) const;

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
         * @param low low region coordinates
         * @param high high region coordinates
         * @param aggs resolver aggregate functions
         * @return region's aggregates, one per aggregate
         */
        IntervalValueVector ComputeAggregate(const Coordinates &low,
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
         * Returns the name of the synopsis.
         *
         * The synopsis name equals the name of the underlying array.
         *
         * @return synopsis name
         */
        const std::string &GetName() const {
            return synopsis_array_->getName();
        }

        /**
         * Preloads this synopsis.
         *
         * Preloading should be performed before concurrent access to the
         * synopsis. If a synopsis is preloaded, threads will not check the
         * validity of the cells and will use them immediately.
         *
         * Since the checking is done via a bool flag, preload must be done
         * at the initialization phase.
         */
        void Preload();

    private:
        /*
         * make RegionIterator a friend since it requires frequent access
         * to meta-data.
         */
        friend class RegionIterator;

        /*
         *  Parses chunk sizes out of the string. The string is supposed to
         *  have the format "x_size,y_size,...".
         */
        void ParseChunkSizes(const std::string &size_param);

        /*
         * Retrieves the synopsis cell from the specified position.
         */
        void FillCellFromArray(const Coordinates &pos, Cell &cell);

        // Init synopsis iterators.
        void InitIterators();

        /*
         *  Retrieve and load the cell iterator is pointing to.
         *    nc_cell is a fallback cell, when cells are not cached
         */
        Cell *GetCurrentCell(const RegionIterator &iter, Cell *nc_cell);

        // Synopsis cells (linearized)
        Cells cells_;

        // Do we cache cells in memory?
        bool cache_cells_ = false;

        // Is this synopses preloaded? No cell validation required, if it is.
        bool preloaded_ = false;

        // Attribute IDs for min/max/count/sum elements in the synopsis array
        AttributeID min_id_, max_id_, count_id_, sum_id_;

        // Iterators to the synopsis array for faster chunk retrieval
        ConstItemIteratorPtr min_it_, max_it_, sum_it_, count_it_;

        // The sample array
        const ArrayPtr synopsis_array_;

        // The size of a cell (one per dimension)
        Coordinates cell_size_;

        // The number of elements in a cell, based on its shape
        uint64_t shape_cell_size_;

        // The starting and ending points of the synopsis (original coordinates)
        Coordinates synopsis_origin_, synopsis_end_;

        // The number of cells per each dimension
        Coordinates cell_nums_;

        // To synchronize access during cell loads
        std::mutex mtx_;
    };

    /**
     * Synopsis pointer.
     */
    using SynopsisPtr = std::unique_ptr<Synopsis>;

    /*
     *  Iterator over the cells of a region. We assume that cells are
     *  laid in the row-major order for the purpose of returning the
     *  linear number of the current cell the iterator points to.
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
        RegionIterator(const Synopsis &synopsis, const Coordinates &low,
                const Coordinates &high) :
                    pos_(low),
                    cell_pos_(-1),
                    region_low_(low),
                    region_high_(high),
                    synopsis_(synopsis),
                    valid_(true) {
            cell_pos_ = GetCellPos();
        }

        // prefix ++
        RegionIterator &operator++() {
            size_t i = pos_.size() - 1;
            while ((pos_[i] += synopsis_.cell_size_[i]) > region_high_[i]) {
                pos_[i] = region_low_[i];
                if (i == 0) {
                    valid_ = false;
                    break;
                }
                i--;
            }

            if (valid_) {
                /*
                 *  Align by the cell boundary. We need this to properly
                 *  detect cell pieces the region covers.
                 */
                pos_[i] -= (pos_[i] - synopsis_.synopsis_origin_[i]) %
                        synopsis_.cell_size_[i];
                if (i == pos_.size() - 1) {
                    cell_pos_++;
                } else {
                    cell_pos_ = GetCellPos();
                }
            }

            return *this;
        }

        // Is the iterator invalid (at the end?)
        bool end() const {
            return !valid_;
        }

        // Return the linear cell position
        uint64_t GetCell() const {
            return cell_pos_;
        }

        // Return current cell position (non-linear)
        const Coordinates &GetCurrentPosition() const {
            return pos_;
        }

        // Are we fully covering the current cell with the region?
        bool CoversFullCell() const {
            bool res = true;
            for (size_t i = 0; i < pos_.size(); i++) {
                const Coordinate chunk_size = synopsis_.cell_size_[i];
                /*
                 * The point must be aligned with the leftmost corner of the
                 * cell + the region must cover it fully.
                 */
                if (pos_[i] % chunk_size != 0 ||
                        region_high_[i] - pos_[i] + 1 < chunk_size) {
                    res = false;
                    break;
                }
            }
            return res;
        }

        /*
         *  Return the covered portion of the current cell. If the
         *  cell is fully covered it will return its coordinates.
         */
        void GetCoveredCell(Coordinates &low, Coordinates &high) const {
            low = pos_;
            high.resize(pos_.size());
            GetHighCellCoords(high);
        }

        /*
         * Computes the number of elements in the covered part of the cell.
         *
         * The number of elements is shape-based, some might be empty if the
         * contents of the cell.
         */
        uint64_t GetPartSize() const {
            Coordinates high_pos(pos_.size());
            GetHighCellCoords(high_pos);
            uint64_t res = 1;
            for (size_t i = 0; i < pos_.size(); i++) {
                res *= (high_pos[i] - pos_[i] + 1);
            }
            return res;
        }

        /**
         * Return current position in synopsis coordinates.
         *
         * Synopsis coordinates means the position is scaled down by
         * using the cell size. These coordinates are suitable to
         * directly access the synopsis array.
         *
         * @return current position in synopsis coordinates
         */
        Coordinates GetCurrentSynopsisPosition() const {
            Coordinates res{pos_};
            for (size_t i = 0; i < res.size(); i++) {
                res[i] = (res[i] - synopsis_.synopsis_origin_[i]) /
                        synopsis_.cell_size_[i];
            }
            return res;
        }

    private:
        // Returns the linear cell number for the current position
        uint64_t GetCellPos() const {
            // Linear position is computed in synopsis coordinates
            Coordinates cell_coord{GetCurrentSynopsisPosition()};

            // TODO: single out one-/two-dimensional cases?
            uint64_t pos = 0;
            for (size_t i = 0; i < cell_coord.size(); i++) {
                pos *= synopsis_.cell_nums_[i];
                pos += cell_coord[i];
            }

            return pos;
        }

        /*
         * Computes rightmost coordinates for the current cell. The array
         * must be of the appropriate size.
         */
        void GetHighCellCoords(Coordinates &high) const {
            for (size_t i = 0; i < pos_.size(); i++) {
                const Coordinate cell_size = synopsis_.cell_size_[i];
                high[i] = (pos_[i] -
                      (pos_[i] - synopsis_.synopsis_origin_[i]) % cell_size) +
                      cell_size - 1;
                if (high[i] > region_high_[i]) {
                    high[i] = region_high_[i];
                }
            }
        }

        // Current position
        Coordinates pos_;

        // Current linear cell position (to use for the cell cache)
        uint64_t cell_pos_;

        // Boundaries for the region
        const Coordinates &region_low_, &region_high_;

        // The synopsis we are traversing
        const Synopsis &synopsis_;

        // Do we point to a valid position?
        bool valid_;
    };

    /*
     * Every synopsis array name consists of:
     *   <user-defined name>_<attr_name>_NxN... where N are synopsis cell
     *   sizes.
     *
     *   This function just returns the three parts in separate strings.
     */
    static StringVector ParseArrayParamsFromName(const std::string &array_name);

    // Synopses catalog (attr. name --> list of synopsis)
    std::unordered_map<std::string,std::vector<ArrayPtr>> array_synopses_;

    /*
     * Catalog of prepared synopses. First vector is indexed by internal
     * attribute ids, second vector contains synopses ranged by the cell
     * size in decreasing order.
     */
    std::vector<std::vector<SynopsisPtr>> synopses_;

    // Descriptor of the data array we store synopses for
    const ArrayDesc data_desc_;

    // Searchlight config
    const SearchlightConfig &sl_config_;

    // Aggregate map to resolve aggregates
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
            uint64_t part_size, const Sampler::Cell &chunk) = 0;

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

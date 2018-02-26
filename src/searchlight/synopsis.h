/*
 * Copyright 2016, Brown University, Providence, RI.
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
 * @file synopsis.h
 * Synopsis basic definitions and declarations.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_SYNOPSIS_H_
#define SEARCHLIGHT_SYNOPSIS_H_

#include "scidb_inc.h"
#include "base.h"

#include <mutex>
#include <boost/tokenizer.hpp>

namespace searchlight {

/**
 * Type of caching for the synopsis.
 */
enum class CachingType {
    EAGER, /** All cells are in a linear array, high footprint. */
    LAZY,  /** All _requested_ chunks are in a set of arrays. */
    NONE   /** No caching. Each cell is read via array iterators.*/
};

/**
 * Grid synopsis basic definition.
 */
template<class Cell, class CellReader>
class GridSynopsis {
public:
    /**
     *  Iterator over the cells of a region.
     *
     *  We iterate over the cells the given region intersects (fully or
     *  partially). We assume that cells are laid out in the row-major order.
     */
    class RegionIterator {
    public:
        /**
         * Creates an iterator over the region, specified by low and
         * high coordinates (both inclusive).
         *
         * We assume that the low-high coordinates comprise a valid region,
         * where high[i] >= low[i].
         */
        RegionIterator(GridSynopsis &synopsis, const Coordinates &low,
                       const Coordinates &high) :
                region_low_(low),
                region_high_(high),
                synopsis_(synopsis),
                cell_reader_{synopsis.synopsis_array_, synopsis.attributes_},
                pos_{low},
                cell_pos_{GetCellPos()} {}

        /**
         * Reset this iterator to a new region.
         *
         * @param low left region coordinates
         * @param high right region coordinates
         */
        void Reset(const Coordinates &low, const Coordinates &high) {
            region_low_ = low;
            region_high_ = high;
            pos_ = low;
            cell_pos_ = GetCellPos();
            valid_ = true;
        }

        /**
         * Prefix ++ operator.
         *
         * @return this iterator
         */
        RegionIterator &operator++() {
            // Row-major order
            size_t i = pos_.size() - 1;
            /*
             * We HAVE to align the coordinate by the cell here. Otherwise,
             * we might miss small pieces of the region protruding beyond the
             * last cell in any dimension.
             */
            while ((pos_[i] = AlignCoordinate(pos_[i], i) +
                    synopsis_.cell_size_[i]) > region_high_[i]) {
                pos_[i] = region_low_[i];
                if (i == 0) {
                    valid_ = false;
                    break;
                }
                i--;
            }

            if (valid_) {
                /*
                 * Small optimization for moving the cell number: since
                 * we're going in the row-major order, we can just increment
                 * it along a row. However, when me move to the next row, we
                 * recompute it from scratch.
                 */
                if (i == pos_.size() - 1) {
                    cell_pos_++;
                } else {
                    cell_pos_ = GetCellPos();
                }
            }

            return *this;
        }

        /**
         * Check if the iterator is finished.
         *
         * @return true, if the iterator is finished; false, otherwise
         */
        bool end() const {
            return !valid_;
        }

        /**
         * Return current cell's linear position.
         *
         * @return current cell's linear position
         */
        uint64_t GetCellLinear() const {
            return cell_pos_;
        }

        /**
         * Return the current cell.
         *
         * @return current cell
         */
        inline const Cell &GetCell() const {
            return synopsis_.GetCell(*this);
        }

        /**
         * Check if the current cell is covered fully by the region.
         *
         * @return true, if covered fully
         */
        bool CoversFullCell() const {
            bool res = true;
            for (size_t i = 0; i < pos_.size(); i++) {
                const Coordinate chunk_size = synopsis_.cell_size_[i];
                /*
                 * The point must be aligned with the leftmost corner of the
                 * cell + the region must cover it fully.
                 */
                const Coordinate pos_offset =
                        pos_[i] - synopsis_.synopsis_origin_[i];
                if (pos_offset % chunk_size != 0 ||
                        region_high_[i] - pos_[i] + 1 < chunk_size) {
                    res = false;
                    break;
                }
            }
            return res;
        }

        /**
         *  Return the covered portion of the current cell.
         *
         *  If the cell is fully covered it will return the entire cell
         *  coordinates.
         *
         *  @param low low point (out)
         *  @param high high point (out)
         */
        void GetCoveredCell(Coordinates &low, Coordinates &high) const {
            low = pos_;
            high.resize(pos_.size());
            GetHighCellCoords(high);
        }

        /**
         * Compute the number of elements in the covered part of the cell.
         *
         * The number of elements is shape-based. If some elements are empty or
         * null, they still count. If the cell if fully covered, the size of
         * the cell will be returned.
         *
         * @return size of the covered portion of the current cell
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
         * Return current position in synopsis array coordinates.
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
        friend class GridSynopsis;

        /**
         * Return cell reader for the iterator.
         *
         * @return cell reader
         */
        CellReader &GetCellReader() const {
            return cell_reader_;
        }

        /*
         * Align the coordinate with the start of the corresponding cell.
         *
         * @param x coordinate to align
         * @param i dimension index
         * @return aligned coordinate
         */
        Coordinate AlignCoordinate(Coordinate x, size_t i) const {
            return x - (x - synopsis_.synopsis_origin_[i]) %
                    synopsis_.cell_size_[i];
        }

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
                high[i] = AlignCoordinate(pos_[i], i) + cell_size - 1;
                if (high[i] > region_high_[i]) {
                    high[i] = region_high_[i];
                }
            }
        }

        // Boundaries for the region
        Coordinates region_low_, region_high_;

        // The synopsis we are traversing
        GridSynopsis &synopsis_;

        // Context to get cells (mutable, it's convenience/cache)
        mutable CellReader cell_reader_;

        // Do we point to a valid position?
        bool valid_ = true;

        // Current position
        Coordinates pos_;

        // Current linear cell position (to use for the cell cache)
        uint64_t cell_pos_;
    };

    /**
     * Create a new synopsis.
     *
     * No cells are loaded at this point. Only meta-data is initialized.
     *
     * Note, the array itself might have multiple dimensions. The user gives
     * the parameter that specifies the number of dimensions to use for the
     * grid itself. The remaining dimensions are considered to be data-specific
     * and should be handled by derived classes.
     *
     * @param data_desc data array descriptor
     * @param array synopsis array
     * @param dims number of grid dimensions
     */
    GridSynopsis(const ArrayDesc &data_desc, const ArrayPtr &array,
        size_t dims) :
            logger_(log4cxx::Logger::getLogger("searchlight.sampler")),
            synopsis_array_(array) {
        // Array descriptor
        const ArrayDesc &synopsis_desc = array->getArrayDesc();

        // By convenience we store sizes in the name after the last '_'
        const std::string &synopsis_config =
               ArrayDesc::makeUnversionedName(synopsis_desc.getName());
        cell_size_.resize(dims);
        const auto params = TokenizeString(synopsis_config, "_");
        if (params.size() < 3) {
            std::ostringstream err_msg;
            err_msg << "Incorrect name for sample array. "
                    "Must be name_attr_NxNx...: name=" << synopsis_config;
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                                   SCIDB_LE_ILLEGAL_OPERATION)
                    << err_msg.str();
        }
        ParseChunkSizes(params.back());

        // Compute sample boundaries (in original coordinates)
        synopsis_origin_.resize(dims);
        synopsis_end_.resize(dims);
        cell_nums_.resize(dims);
        chunk_nums_.resize(dims);
        chunk_sizes_.resize(dims);
        for (size_t i = 0; i < dims; i++) {
            // Metadata is filled from the data descriptor
            const DimensionDesc &data_dim = data_desc.getDimensions()[i];
            /*
             * Cannot use low/high boundary, since SciDB aligns them by the
             * chunk, so high_boundary might not be equal to current_end.
             *
             * See trim() in ArrayDesc for unbounded arrays.
             *
             * Use StartMin() here since it determines the total left boundary,
             * which we need for proper origin.
             */
            synopsis_origin_[i] = data_dim.getStartMin();
            synopsis_end_[i] = data_dim.getCurrEnd();

            // Number of cells (cell sizes are read from the array's name)
            const uint64_t curr_length = synopsis_end_[i] -
                    synopsis_origin_[i] + 1;
            cell_nums_[i] = (curr_length - 1) / cell_size_[i] + 1;

            // Synopsis correctness checking
            const DimensionDesc &syn_dim = synopsis_desc.getDimensions()[i];
            if (syn_dim.getStartMin() != 0) {
                throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                        SCIDB_LE_ILLEGAL_OPERATION)
                        << "Synopsis coordinates must start from 0";
            }
            // Synopsis chunks info
            chunk_sizes_[i] = syn_dim.getChunkInterval();
            chunk_nums_[i] = syn_dim.getCurrEnd() / chunk_sizes_[i] + 1;
        }

        /*
         * Find min/max/sum/count ids.
         */
        for (const auto &attr: synopsis_desc.getAttributes(true)) {
            attributes_.emplace(attr.getName(), attr.getId());
        }
    }

    /**
     * Destructor.
     */
    virtual ~GridSynopsis() = default;

    /**
     * Set caching type for the synopsis.
     *
     * @param mode see CachingType
     */
    void SetCacheType(CachingType mode) {
        cache_type_ = mode;
        if (mode == CachingType::EAGER) {
            /*
             * Create the entire cell cache right now to simplify concurrency
             * later. The cells are invalid and will be loaded from the synopsis
             * later in lazy fashion.
             */
            chunks_ = Chunks(1);
            chunks_[0].cells_ = Cells(GetTotalCellCount());
            chunks_[0].valid_.store(true, std::memory_order_release);
        } else if (mode == CachingType::LAZY) {
            /*
             * Reserve space for array chunks. All chunks are invalid at the
             * beginning.
             */
            chunks_ = Chunks(GetTotalChunksCount());
        } else {
            chunks_ = Chunks(); // make it empty
        }
    }

    /**
     * Check if this synopsis is cached.
     *
     * @return true, if the cells are cached; false, otherwise
     */
    bool IsCached() const {
        return cache_type_ != CachingType::NONE;
    }

    /**
     * Return cahing type for the synopsis.
     *
     * @return synopsis caching type
     */
    CachingType GetCachingType() const {
        return cache_type_;
    }

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
    void CheckAndCorrectBounds(Coordinates &low, Coordinates &high) const {
        const size_t dims = synopsis_origin_.size();
        for (size_t i = 0; i < dims; i++) {
            if (low[i] < synopsis_origin_[i]) {
                low[i] = synopsis_origin_[i];
            }
            if (high[i] > synopsis_end_[i]) {
                high[i] = synopsis_end_[i];
            }
        }
    }

    /**
     * Check if the region is valid for the synopsis.
     *
     * The region is valid if it intersects with the synopsis region.
     *
     * @param low low coordinates of the region
     * @param high high coordinates of the region
     * @return true, if valid; false, otherwise
     */
    bool CheckIfValid(const Coordinates &low, const Coordinates &high) const {
        const size_t dims = synopsis_origin_.size();
        if (low.size() != high.size() || low.size() != dims) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                    SCIDB_LE_ILLEGAL_OPERATION)
                    << "Specified region has inconsistent dimensions!";
        }
        for (size_t i = 0; i < dims; ++i) {
            if (low[i] > high[i]) {
                throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                        SCIDB_LE_ILLEGAL_OPERATION)
                        << "Specified region has inconsistent dimensions!";
            }
            if (high[i] < synopsis_origin_[i] || low[i] > synopsis_end_[i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * Check if the point's coordinates are within the synopsis area.
     *
     * @param point point's coordinates
     * @return true, if the point is within the area; false, otherwise
     */
    bool CheckBounds(const Coordinates &point) const {
        if (point.size() != synopsis_origin_.size()) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                                   SCIDB_LE_ILLEGAL_OPERATION)
                    << "Specified point has inconsistent dimensions!";
        }

        for (size_t i = 0; i < synopsis_origin_.size(); i++) {
            if (point[i] < synopsis_origin_[i] || point[i] > synopsis_end_[i]) {
                return false;
            }
        }

        return true;
    }

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
    void Preload() {
        if (cache_type_ != CachingType::EAGER) {
            LOG4CXX_WARN(logger_, "Attempting to preload non-eager synopsis.");
            return;
        }
        LOG4CXX_INFO(logger_, "Preloading synopsis: " << GetName());
        // preload by touching every cell
        RegionIterator iter{*this, synopsis_origin_, synopsis_end_};
        while (!iter.end()) {
            iter.GetCell();
            ++iter;
        }
        preloaded_ = true;
    }

    /**
     * Return memory size needed to cache all cells of this synopsis.
     *
     * @return required cache memory size (in bytes)
     */
    size_t MemorySize() const {
        return MemoryPerCell() * GetTotalCellCount();
    }

    /**
     * Output this synopsis statistics into a stream.
     *
     * @param str stream for output
     */
    void OutputStats(std::ostream &str) const {
        str << "Synopsis " << GetName() << ":\n";
        str << "\tCells accessed: "
                << cells_accessed_.load(std::memory_order_relaxed);
        str << "\n\tMemory footprint (MB): "
                << ComputeMemoryFootprint() / 1024 / 1024;
    }

    /**
     * Return MBR for the region aligned with the synopsis grid.
     *
     * MBR for the region means the following: the original region is extended
     * until it intersects all cells completely (there are no partial overlaps).
     * MBR is the minimal extension.
     *
     * @param low  low region boundary (in/out)
     * @param high high region boundary (in/out)
     */
    void GetSynopsisMBR(Coordinates &low, Coordinates &high) const {
        for (size_t i = 0; i < low.size(); ++i) {
            low[i] -= (low[i] - synopsis_origin_[i]) % cell_size_[i];
            high[i] -= (high[i] - synopsis_origin_[i]) % cell_size_[i];
            high[i] += cell_size_[i] - 1;
        }
    }

    /**
     * Return the cost of a region in cells.
     *
     * @param reg the region to compute the cost for
     * @return region cost
     */
    size_t GetRegionCost(const Coordinates &low,
            const Coordinates &high) const {
        /*
         * We align the region with the grid to get proper estimations. Note,
         * it might be possible to go without it, by counting the even
         * number of cells and adding 1 for the remainder, if any. However,
         * such an estimation might be off by 1 for every dimension,
         * if the remainders are on both ends of the interval.
         */
        Coordinates low_mbr{low}, high_mbr{high};
        GetSynopsisMBR(low_mbr, high_mbr);
        size_t cost = 1;
        for (size_t i = 0; i < low_mbr.size(); ++i) {
            const size_t len = high_mbr[i] - low_mbr[i] + 1;
            assert(len % cell_size_[i] == 0);
            cost *= len / cell_size_[i];
        }
        return cost;
    }

    /**
     * Replace physical synopsis array with the new one.
     *
     * This is actually useful for synopses cached in the inter-query
     * cache, since arrays are tied to scidb::Query.
     *
     * Note, the cell cache remains intact.
     *
     * @param array new synopsis array
     */
    void ReplaceArray(const ArrayPtr &array) {
        synopsis_array_ = array;
    }

private:
    /*
     * Memory required to store each cell.
     *
     * Cells might have some additional info attached, so the function
     * can be overriden to reflect that.
     *
     * @return memory per cell
     */
    virtual size_t MemoryPerCell() const {
        return sizeof(Cell);
    }

    // Vector of cells
    using Cells = std::vector<Cell>;

    // Synopsis chunk
    struct SynopsisChunk {
        std::atomic<bool> valid_;
        Cells cells_;

        // Default validity constructor
        SynopsisChunk(bool valid = false) : valid_{valid} {}

        // Copy constructor (needed to copy atomic)
        SynopsisChunk(SynopsisChunk &&other) :
            valid_{other.valid_.load(std::memory_order_relaxed)},
            cells_{std::move(other.cells_)} {}
    };
    // Synopsis chunks
    using Chunks = std::vector<SynopsisChunk>;

    /*
     * Structure to hold synopsis cell coordinates, where the synopsis is
     * divided into chunks. Thus, cell coordinates, both vector and linear,
     * are respective to the chunk's origin.
     */
    struct ChunkCellCoordinates {
        Coordinates chunk_;
        size_t chunk_linear_ = 0;
        Coordinates cell_;
        size_t cell_chunk_linear_ = 0;

        /*
         * Init coordinates from cell coordinates.
         *
         * Synposis needed to compute other coordinates
         */
        ChunkCellCoordinates(const Coordinates &cell_coords,
                const GridSynopsis &syn) :
                    chunk_(cell_coords.size()),
                    cell_{cell_coords} {
            for (size_t i = 0; i < cell_.size(); ++i) {
                // Chunk coords
                chunk_[i] = cell_[i] - cell_[i] % syn.chunk_sizes_[i];
                chunk_linear_ *= syn.chunk_nums_[i];
                chunk_linear_ += chunk_[i] / syn.chunk_sizes_[i];
                // Cell coords
                cell_chunk_linear_ *= syn.chunk_sizes_[i];
                cell_chunk_linear_ += cell_[i] - chunk_[i];
            }
        }
    };

    /**
     * Return cell pointed by the iterator.
     *
     * The user specifies contex, containing additional information
     * for the reader.
     *
     * @param iter iterator pointing at the cell
     * @param ctx reader context
     * @return cell pointed to by the iterator
     */
    const Cell &GetCell(const RegionIterator &iter) {
        CellReader &reader = iter.GetCellReader();
        if (cache_type_ == CachingType::EAGER) {
            assert(chunks_.size() == 1 &&
                    iter.GetCellLinear() < chunks_[0].cells_.size());
            Cell &cell = chunks_[0].cells_[iter.GetCellLinear()];

            // Load cell if needed (cannot check validity here -- concurrency!)
            if (!preloaded_) {
                /*
                 * Need a lock here. For efficiency reasons this has been
                 * implemented as a double-checked locking with fences and atomics.
                 *
                 * This used to be a separate function, but later it was inlined
                 * to avoid a huge overhead on non-preloaded synopses. This is a
                 * critical code path!
                 */
                bool valid = cell.valid_.load(std::memory_order_acquire);
                if (!valid) {
                    std::lock_guard<std::mutex> lock{mtx_};
                    valid = cell.valid_.load(std::memory_order_relaxed);
                    if (!valid) {
                        reader.Read(iter.GetCurrentSynopsisPosition(), cell);
                        cell.valid_.store(true, std::memory_order_release);
                    }
                }
            }
            assert(cell.valid_);
            return cell;
        } else if (cache_type_ == CachingType::LAZY) {
            // Determine position
            ChunkCellCoordinates pos{iter.GetCurrentSynopsisPosition(), *this};
            SynopsisChunk &chunk = chunks_[pos.chunk_linear_];

            // Fetch (see the double-checked locking comment above)
            bool valid = chunk.valid_.load(std::memory_order_acquire);
            if (!valid) {
                std::lock_guard<std::mutex> lock{mtx_};
                valid = chunk.valid_.load(std::memory_order_relaxed);
                if (!valid) {
                    LOG4CXX_TRACE(logger_, "Lazy-loading chunk " <<
                                  pos.chunk_linear_);
                    FillCellsFromChunk(pos, reader);
                    chunk.valid_.store(true, std::memory_order_release);
                }
            }

            // The cell
            assert(chunk.cells_[pos.cell_chunk_linear_].valid_);
            return chunk.cells_[pos.cell_chunk_linear_];
        } else {
            reader.Read(iter.GetCurrentSynopsisPosition(), reader.TempCell());
            return reader.TempCell();
        }
    }

    /*
     *  Parses chunk sizes out of the string. The string is supposed to
     *  have the format "x_size,y_size,...".
     */
    void ParseChunkSizes(const std::string &size_param) {
        using TokenSeparator = boost::char_separator<char>;
        using Tokenizer = boost::tokenizer<TokenSeparator>;
        TokenSeparator sep("xX"); // size_1xsize_2x...xsize_n
        Tokenizer tokenizer(size_param, sep);

        shape_cell_size_ = 1;
        size_t i = 0;
        for (auto cit = tokenizer.begin(); cit != tokenizer.end(); ++cit) {
            cell_size_[i] = std::stoi(*cit);
            assert(cell_size_[i] > 0);
            shape_cell_size_ *= cell_size_[i];
            i++;
        }

        if (i != cell_size_.size()) {
            std::ostringstream err_msg;
            err_msg << "Could not retrieve all cell sizes: conf=" <<
                    size_param << ", needed sizes=" << cell_size_.size();
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << err_msg.str();
        }
    }

    /**
     * Return total number of cells in the synopsis.
     *
     * @return total number of synopsis cells
     */
    size_t GetTotalCellCount() const {
        size_t total_cell_count = 1;
        for (auto cn: cell_nums_) {
            total_cell_count *= cn;
        }
        return total_cell_count;
    }

    /**
     * Return total number of synopsis array chunks.
     *
     * @return total number of synopsis array chunks
     */
    size_t GetTotalChunksCount() const {
        size_t total_chunks_count = 1;
        for (auto cn: chunk_nums_) {
            total_chunks_count *= cn;
        }
        return total_chunks_count;
    }

    /*
     * Fill all cells from the specified chunk.
     */
    void FillCellsFromChunk(const ChunkCellCoordinates &pos,
            CellReader &reader) const {
        const size_t chunk_index = pos.chunk_linear_;
        assert(chunk_index < chunks_.size() && !chunks_[chunk_index].valid_);
        SynopsisChunk &chunk = chunks_[chunk_index];

        // The number of cells in a chunk
        size_t total_cells_per_chunk = 1;
        for (auto cs: chunk_sizes_) {
            total_cells_per_chunk *= cs;
        }
        // Will make it valid later (see GetCell()).
        chunk.cells_ = Cells(total_cells_per_chunk);

        // Start and end chunk positions
        const Coordinates &chunk_start = pos.chunk_;
        Coordinates chunk_end(chunk_start);
        for (size_t i = 0; i < chunk_start.size(); ++i) {
            chunk_end[i] += chunk_sizes_[i] - 1;
        }

        // Retrieve the chunk's cells
        Coordinates current_pos{pos.chunk_};
        bool cell_is_valid = true;
        while (cell_is_valid) {
            // Get the cell and make it valid
            ChunkCellCoordinates chunk_cell_pos{current_pos, *this};
            Cell &cell = chunk.cells_[chunk_cell_pos.cell_chunk_linear_];
            reader.Read(chunk_cell_pos.cell_, cell);
            // Make it valid. Relaxed is fine, everything is locked for LAZY.
            cell.valid_.store(true, std::memory_order_relaxed);

            // increment the position
            size_t i = current_pos.size() - 1;
            while (++current_pos[i] > chunk_end[i]) {
                current_pos[i] = chunk_start[i];
                if (i == 0) {
                    cell_is_valid = false;
                    break;
                }
                --i;
            }
        }
    }

    // Retun memory footprint for the synopsis
    size_t ComputeMemoryFootprint() const {
        if (cache_type_ == CachingType::EAGER) {
            return MemorySize() + sizeof(SynopsisChunk);
        } else if (cache_type_ == CachingType::LAZY) {
            size_t res = GetTotalChunksCount() * sizeof(SynopsisChunk);
            /*
             *  Count valid chunks' cells. For some cells the estimation
             *  may be quite off, if cells have some additional info beyond
             *  there sizeof(). This is because some synopsis cells might be
             *  just empty.
             */
            const size_t cell_mem = MemoryPerCell();
            for (const auto &c: chunks_) {
                if (c.valid_.load(std::memory_order_relaxed)) {
                    res += c.cells_.size() * cell_mem;
                }
            }
            return res;
        } else {
            // None
            return 0;
        }
    }

    /*
     *  Synopsis cache. Usage depends on the caching type:
     *
     *  EAGER: there is one entry [0] for all synopsis cells.
     *  LAZY: entries correspond to synopsis array chunks.
     */
    mutable Chunks chunks_;  // Synopsis chunks with linearized cells

    // Type of cache for the synopsis
    CachingType cache_type_ = CachingType::NONE;

    // Is this synopses preloaded? No cell validation required, if it is.
    bool preloaded_ = false;

    /* To synchronize access during cell loads.
     *
     * TODO: Strictly speaking, there is no need to have a single mutex for
     * multiple consumers loading different cells. Another simple solution
     * is to have two mutexes, one for even and another for odd cells.
     * One could even generalize further, e.g., to four mutexes. This way
     * if consumers load different (modulo) cells, they won't block and
     * if they load the same cell, one of them will block, which is
     * the main intention here. Note, in this case each consumer should
     * have its own copy of the iterators.
     */
    mutable std::mutex mtx_;

    /*
     * Need a copy of the logger here. We use the one from the sampler.
     */
    log4cxx::LoggerPtr logger_;

protected:
    // Stats: total number of cell accessed
    std::atomic<uint64_t> cells_accessed_{0};

    // The sample array
    ArrayPtr synopsis_array_;

    // The size of a cell (one per dimension)
    Coordinates cell_size_;

    // The number of elements in a cell, based on its shape
    uint64_t shape_cell_size_;

    // The starting and ending points of the synopsis (original coordinates)
    Coordinates synopsis_origin_, synopsis_end_;

    // The number of cells per each dimension
    Coordinates cell_nums_;

    // Number of synopsis array chunks in each dimension
    SizeVector chunk_nums_;

    // Size of the synopsis array chunks for each dimension
    SizeVector chunk_sizes_;

    // Attribute map
    AttributeMap attributes_;
};

/**
 * A region (cell) of the synopsis.
 */
struct AggCell {
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
    AggCell(double min, double max, double sum, uint64_t count) :
        min_(min), max_(max), sum_(sum), count_(count), valid_{true} {}

    /**
     * Copy constructor.
     *
     * Needed because of the atomic and usage of this class in a vector.
     *
     * @param cell cell to copy from
     */
    AggCell(const AggCell &cell) :
        min_(cell.min_), max_(cell.max_), sum_(cell.sum_),
        count_(cell.count_),
        valid_(cell.valid_.load(std::memory_order_relaxed)) {}

    /**
     * Constructs an invalid (default) cell.
     *
     * If the cell is invalid, that means it has to be loaded from disk.
     */
    AggCell() = default;
};

/**
 * Sequence synopsis cell.
 *
 * Basically trace MBRs with low/high coordinates.
 */
struct SeqCell {
    /**
     * MBR for storing a trace region.
     */
    struct {
        DoubleVector low_, high_;

        /**
         * Compute minimum square distance from point to this rectangle.
         *
         * @param point point coordinates
         * @return minimum square distance
         */
        double MinSqDist(const double *point) const {
            assert(low_.size() == high_.size());

            double res = 0;
            for (size_t i = 0; i < low_.size(); ++i) {
                double min_coord_dist = std::max(point[i] - high_[i], 0.0);
                min_coord_dist = std::max(min_coord_dist, low_[i] - point[i]);
                res += min_coord_dist * min_coord_dist;
            }
            return res;
        }
    } mbr_;

    /**
     * Cell has been loaded.
     */
    std::atomic<bool> valid_{false};

    /**
     * Copy constructor.
     *
     * Needed because of the atomic and usage of this class in a vector.
     *
     * @param cell cell to copy from
     */
    SeqCell(SeqCell &&cell) :
        mbr_(std::move(cell.mbr_)),
        valid_(cell.valid_.load(std::memory_order_relaxed)) {}

    /**
     * Constructs an invalid (default) cell.
     *
     * If the cell is invalid, that means it has to be loaded from disk.
     */
    SeqCell() = default;
};

/**
 * Aggregate cell reader for synopsis arrays.
 */
class AggCellItemReader {
public:
    // Constructor.
    AggCellItemReader(const ArrayPtr &array, const AttributeMap &attrs) :
        array_{array},
        attributes_(attrs) {}

    // Reads a cell
    void Read(const Coordinates &pos, AggCell &cell);

    // Temporary cell
    AggCell &TempCell() {
        return cell_;
    }

private:
    // Init iterators
    void InitIterators();

    ArrayPtr array_;
    const AttributeMap &attributes_;
    ConstItemIteratorPtr min_it_, max_it_, sum_it_, count_it_;
    AggCell cell_;
};

/**
 * Seqcell reader for synopsis arrays.
 */
class SeqCellItemReader {
public:
    // Constructor.
    SeqCellItemReader(const ArrayPtr &array, const AttributeMap &attrs) :
        array_{array},
        attributes_(attrs) {}

    // Reads a cell
    void Read(const Coordinates &pos, SeqCell &cell);

    // Temporary cell
    SeqCell &TempCell() {
        return cell_;
    }

private:
    // Init iterators
    void InitIterators();

    ArrayPtr array_;
    const AttributeMap &attributes_;
    ConstItemIteratorPtr low_it_, high_it_;
    size_t coords_num_ = 0; // number of DFT components
    Coordinates read_pos_;
    SeqCell cell_;
};

/**
 * Info about the transformed sequence.
 */
struct TransformedSequenceInfo {
    // Original sequence length
    size_t original_length_;
    // Attribute of the sequence (the search one)
    AttributeID sattr_;
    // The transformed sequences (subsequence size --> transformed)
    std::unordered_map<size_t, DoubleVector> sequence_;

    TransformedSequenceInfo(size_t len, AttributeID sattr) :
        original_length_(len), sattr_(sattr) {}
};
} /* namespace searchlight */
#endif /* SEARCHLIGHT_SYNOPSIS_H_ */

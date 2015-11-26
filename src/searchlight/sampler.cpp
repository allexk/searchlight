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
 * @file sampler.cpp
 * The implementation of the sampler.
 *
 * @author Alexander Kalinin
 */

#include "sampler.h"
#include "array_desc.h"

#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>
#include <boost/tokenizer.hpp>

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.sampler"));

/*
 * TODO: approximate values are now computed... incorrectly. They might be
 * plausible, but a better way would be to compute the probability and
 * possibly return NULL, if suitable. For now, the value will always be
 * something non-NULL, but the error might be quite significant.
 */
class AverageSampleAggregate : public SampleAggregate {
public:
    static SampleAggregate *Create() {
        return new AverageSampleAggregate;
    }

    virtual void AccumulateChunk(uint64_t cell_size, uint64_t part_size,
            const Sampler::Cell &cell) {
        // chunk info
        const double sum = cell.sum_;
        const uint64_t count = cell.count_;
        const double min = cell.min_;
        const double max = cell.max_;

        // is it a full chunk? guaranteed to be non-empty
        if (part_size == cell_size) {
            not_null_ = true;
            sum_ += sum;
            approx_sum_ += sum;
            count_ += count;
            approx_count_ += count;
            return;
        }

        // a part of a chunk; guaranteed non-empty
        const size_t cell_num = cell_min_count_.size();
        const uint64_t empty_count = cell_size - count;

        // compute counts
        const uint64_t min_count = part_size >= empty_count ?
                part_size - empty_count : 0;
        const uint64_t max_count = part_size >= count ? count : part_size;
        if (min_count > 0) {
            not_null_ = true;
        }
        cell_min_count_.push_back(min_count);
        cell_max_count_.push_back(max_count);

        /*
         * For the approximate computation we make the uniformity assumption,
         * including the uniform spread of empty elements. In practice,
         * it is probably less probable -- empty elements tend to stick
         * together. But without any additional knowledge we cannot do better.
         */
        const double part_ratio = double(part_size) / cell_size;
        approx_sum_ += sum * part_ratio;
        approx_count_ += round(part_ratio * count);

        /*
         *  Now, we will figure out the best distributions for upper/lower
         *  boundaries for avg().
         */
        if (max != min) {
            // upper boundary efforts
            uint64_t count_left = count;
            const uint64_t max_num = floor((sum - count * min) / (max - min));
            cell_efforts_upper_.push_back(Effort(max, max_num, cell_num));
            count_left -= max_num;
            const double deltau = sum - max_num * max -
                    (count_left - 1) * min;
            cell_efforts_upper_.push_back(Effort(deltau, 1, cell_num));
            count_left--;
            if (count_left > 0) {
                cell_efforts_upper_.push_back(Effort(min, count_left,
                        cell_num));
            }

            // lower boundary efforts
            count_left = count;
            const uint64_t min_num = floor((count * max - sum) / (max - min));
            cell_efforts_lower_.push_back(Effort(min, min_num, cell_num));
            count_left -= min_num;
            const double deltal = sum - min_num * min -
                    (count_left - 1) * max;
            cell_efforts_lower_.push_back(Effort(deltal, 1, cell_num));
            count_left--;
            if (count_left > 0) {
                cell_efforts_lower_.push_back(Effort(max, count_left,
                        cell_num));
            }
        } else {
            cell_efforts_upper_.push_back(Effort(max, count, cell_num));
            cell_efforts_lower_.push_back(Effort(min, count, cell_num));
        }
    }

    virtual void Finalize(IntervalValue &res) {
        if (not_null_) {
            res.state_ = IntervalValue::NON_NULL;
        } else {
            if (!cell_min_count_.empty()) {
                // got possible non-empty parts
                res.state_ = IntervalValue::MAY_NULL;
            } else {
                res.state_ = IntervalValue::NUL;
                return;
            }
        }

        // we have to make copies of counts to reuse them for the lower bound
        std::vector<uint64_t> chunk_min_count_copy(cell_min_count_);
        std::vector<uint64_t> chunk_max_count_copy(cell_max_count_);

        // sort the efforts by the value
        std::sort(cell_efforts_upper_.begin(), cell_efforts_upper_.end(),
                Effort::GreaterOrder()); // descending
        std::sort(cell_efforts_lower_.begin(), cell_efforts_lower_.end());

        /*
         * Compute the upper bound: raise the average while you can, and then
         * lower as less as possible.
         */
        double current_sum = sum_;
        double current_count = count_;
        for (size_t i = 0; i < cell_efforts_upper_.size(); i++) {
            const Effort &eff = cell_efforts_upper_[i];
            uint64_t &min_count_chunk = cell_min_count_[eff.cell_];
            uint64_t &max_count_chunk = cell_max_count_[eff.cell_];
            if (current_count == 0) {
                FitNumAndSubtractChunkCounts(min_count_chunk,
                    max_count_chunk, eff.number_, current_sum,
                    current_count, eff.value_, true);
            } else {
                if (eff.value_ >= current_sum / current_count) {
                    if (max_count_chunk > 0) {
                        FitNumAndSubtractChunkCounts(min_count_chunk,
                            max_count_chunk, eff.number_, current_sum,
                            current_count, eff.value_, true);
                    }
                } else { // lowering the avg as less as possible
                    if (min_count_chunk > 0) {
                        FitNumAndSubtractChunkCounts(min_count_chunk,
                            max_count_chunk, eff.number_, current_sum,
                            current_count, eff.value_, false);
                    }
                }
            }
        }
        res.max_ = current_sum / current_count; /* current_count > 0 here */

        /*
         * Computer the lower bound: lower the average while you can and
         * then raise it as little as possible.
         */
        current_sum = sum_;
        current_count = count_;
        for (size_t i = 0; i < cell_efforts_lower_.size(); i++) {
            const Effort &eff = cell_efforts_lower_[i];
            uint64_t &min_count_chunk = chunk_min_count_copy[eff.cell_];
            uint64_t &max_count_chunk = chunk_max_count_copy[eff.cell_];
            if (current_count == 0) {
                FitNumAndSubtractChunkCounts(min_count_chunk,
                    max_count_chunk, eff.number_, current_sum,
                    current_count, eff.value_, true);
            } else {
                if (eff.value_ <= current_sum / current_count) {
                    if (max_count_chunk > 0) {
                        FitNumAndSubtractChunkCounts(min_count_chunk,
                            max_count_chunk, eff.number_, current_sum,
                            current_count, eff.value_, true);
                    }
                } else { // raising the avg as less as possible
                    if (min_count_chunk > 0) {
                        FitNumAndSubtractChunkCounts(min_count_chunk,
                            max_count_chunk, eff.number_, current_sum,
                            current_count, eff.value_, false);
                    }
                }
            }
        }
        res.min_ = current_sum / current_count; /* current_count > 0 here */

        // approximate value
        if (approx_count_ > 0) {
            res.val_ = approx_sum_ / approx_count_;
        }
        if (approx_count_ == 0 || res.val_ < res.min_) {
            res.val_ = res.min_;
        }
    }

private:
    AverageSampleAggregate() :
        not_null_(false), sum_(0), count_(0),
        approx_sum_(0), approx_count_(0) {}

    /*
     * Tries to fit num elements inside the count interval and returns
     * the number of elements p really possible to fit (p<=num). base_max
     * means=true means we have as many as max_count elements available,
     * false -- min_count. The counters are modified appropriately to fit the
     * resulting number. sum and count are modified as well by using the
     * p elem values.
     */
    static void FitNumAndSubtractChunkCounts(uint64_t &min_count,
            uint64_t &max_count, uint64_t num, double &sum, double &count,
            double elem, bool base_max) {
        uint64_t fit = base_max ? std::min(num, max_count) :
                std::min(num, min_count);

        max_count -= fit;
        min_count = min_count >= fit ? min_count - fit : 0;

        sum += fit * elem;
        count += fit;
    }

    /*
     *  An effort describes the cell's best effort to bring up/down the avg.
     *  For example, to increase the avg. maximally we would add the value_
     *  the number_ times to the set of values composing the upper bound.
     *  Cell_ references the cell number for this effort.
     */
    struct Effort {
        double value_;
        uint64_t number_;
        size_t cell_;

        Effort(double value, uint64_t number, size_t cell) :
            value_(value), number_(number), cell_(cell) {}

        bool operator<(const Effort &other) const {
            return value_ < other.value_;
        }

        struct GreaterOrder {
            bool operator()(const Effort &e1, const Effort &e2) const {
                return e1.value_ > e2.value_;
            }
        };
    };

    // Best efforts for the upper and lower bounds
    std::vector<Effort> cell_efforts_upper_, cell_efforts_lower_;

    /*
     * The minimum and maximum number of elements for each cell.
     */
    std::vector<uint64_t> cell_min_count_, cell_max_count_;

    // Is it definitely not null?
    bool not_null_;

    // Current sum and count for boundaries computation
    double sum_;
    uint64_t count_;

    // sum and count for the approximate value
    double approx_sum_;
    uint64_t approx_count_;
};

class SumSampleAggregate : public SampleAggregate {
public:
    static SampleAggregate *Create() {
        return new SumSampleAggregate;
    }

    virtual void AccumulateChunk(uint64_t cell_size,
            uint64_t part_size, const Sampler::Cell &cell) {
        // cell info
        const double sum = cell.sum_;
        const uint64_t count = cell.count_;
        const double min = cell.min_;
        const double max = cell.max_;

        // it cannot be null, since cells are non-empty here
        null_ = false;

        // is it a full cell? guaranteed to be non-empty
        if (part_size == cell_size) {
            not_null_ = true;
            min_sum_ += sum;
            max_sum_ += sum;
            approx_sum_ += sum;
            return;
        }

        // compute counts
        const uint64_t empty_count = cell_size - count;
        const uint64_t min_count = part_size >= empty_count ?
                part_size - empty_count : 0;
        const uint64_t max_count = part_size >= count ? count : part_size;
        if (min_count > 0) {
            not_null_ = true;
        }

        // trying to get better bounds by using sum, instead of trivial ones
        const double min_sum = std::max(sum - (count - min_count) * max,
                min_count * min);
        const double max_sum = std::min(sum - (count - max_count) * min,
                max * max_count);
        min_sum_ += min_sum;
        max_sum_ += max_sum;

        /*
         * For the approximate computation we make the uniformity assumption,
         * including the uniform spread of empty elements. In practice,
         * it is probably less probable -- empty elements tend to stick
         * together. But without any additional knowledge we cannot do better.
         */
        const double part_ratio = double(part_size) / cell_size;
        approx_sum_ += sum * part_ratio;
    }

    virtual void Finalize(IntervalValue &res) {
        if (not_null_) {
            res.state_ = IntervalValue::NON_NULL;
        } else {
            if (!null_) {
                // got possible non-empty parts
                res.state_ = IntervalValue::MAY_NULL;
            } else {
                res.state_ = IntervalValue::NUL;
                return;
            }
        }

        res.min_ = min_sum_;
        res.max_ = max_sum_;
        res.val_ = approx_sum_ < min_sum_ ? min_sum_ : approx_sum_;
    }

private:
    SumSampleAggregate() :
        min_sum_(0), max_sum_(0), approx_sum_(0),
        not_null_(false), null_(true) {}

    // Minimum/maximum and approximate sum
    double min_sum_, max_sum_, approx_sum_;

    // is it definitely not null? or is it definitely null?
    bool not_null_, null_;
};

class MinMaxSampleAggregate : public SampleAggregate {
public:
    static SampleAggregate *CreateMin() {
        return new MinMaxSampleAggregate(true);
    }

    static SampleAggregate *CreateMax() {
        return new MinMaxSampleAggregate(false);
    }

    virtual void AccumulateChunk(uint64_t cell_size,
            uint64_t part_size, const Sampler::Cell &cell) {
        // cell info
        const double min = cell.min_;
        const double max = cell.max_;

        // it cannot be null, since cells are non-empty here
        null_ = false;

        // is it a full cell? guaranteed to be non-empty
        if (part_size == cell_size) {
            not_null_ = true;
            // here we know the exact min and max
            if (is_min_) {
                max_ = std::min(max_, min);
                min_ = std::min(min_, min);
            } else {
                max_ = std::max(max_, max);
                min_ = std::max(min_, max);
            }
            return;
        }

        // compute counts
        const uint64_t empty_count = cell_size - cell.count_;
        const uint64_t min_count = part_size >= empty_count ?
                part_size - empty_count : 0;
        if (min_count > 0) {
            not_null_ = true;
        }

        // partial match: a range of values
        if (is_min_) {
            min_ = std::min(min_, min);
            max_ = std::min(max_, max);
        } else {
            min_ = std::max(min_, min);
            max_ = std::max(max_, max);
        }
    }

    virtual void Finalize(IntervalValue &res) {
        if (not_null_) {
            res.state_ = IntervalValue::NON_NULL;
        } else {
            if (!null_) {
                // got possible non-empty parts
                res.state_ = IntervalValue::MAY_NULL;
            } else {
                res.state_ = IntervalValue::NUL;
                return;
            }
        }

        res.min_ = min_;
        res.max_ = max_;
        /*
         * Approximate value: the middle point, which guarantees the best
         * worst-case error.
         */
        res.val_ = min_ + (max_ - min_) / 2;
    }

private:
    MinMaxSampleAggregate(bool min) :
        approx_(0),
        is_min_(min), not_null_(false), null_(true) {
        if (is_min_) {
            min_ = max_ = std::numeric_limits<double>::max();
        } else {
            min_ = max_ = std::numeric_limits<double>::lowest();
        }
    }

    // Minimum/maximum and approximate values
    double min_, max_, approx_;

    // is it the min aggregate?
    bool is_min_;

    // is it definitely not null? or is it definitely null?
    bool not_null_, null_;
};

Sampler::DFTSynopsis::DFTSynopsis(const ArrayDesc &data_desc,
        const ArrayPtr &array) : synopsis_array_(array) {
    // Array descriptor
    const ArrayDesc &synopsis_desc = array->getArrayDesc();
    // By convenience we store sizes in the name after the last '_'
    const std::string &synopsis_config =
            ArrayDesc::makeUnversionedName(synopsis_desc.getName());
    ParseDFTSynopsisParams(ParseArrayParamsFromName(synopsis_config).back());

    // First dimension is trace MBRs: parse synopsis shape from there
    const DimensionDesc &data_dim = data_desc.getDimensions()[0];
    /*
     * Cannot use low/high boundary, since SciDB aligns them by the
     * chunk, so high_boundary might not be equal to current_end.
     *
     * See trim() in ArrayDesc for unbounded arrays.
     *
     * Use StartMin() here since it determines the total left boundary,
     * which we need for proper origin.
     */
    synopsis_origin_ = data_dim.getStartMin();
    // Covered sequences must fit exactly, should be careful with synopsis_end_
    synopsis_end_ = data_dim.getCurrEnd() - cell_size_ + 1;
    // Number of cells (cell sizes are read from the array's name)
    const uint64_t curr_length = synopsis_end_ - synopsis_origin_ + 1;
    cell_num_ = (curr_length - 1) / cell_size_ + 1;
    // Synopsis correctness checking
    const DimensionDesc &syn_dim = synopsis_desc.getDimensions()[0];
    if (syn_dim.getStartMin() != 0) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Synopsis coordinates must start from 0";
    }
    if (syn_dim.getCurrEnd() + 1 < cell_num_) { // StartMin() == 0
        std::ostringstream err_msg;
        err_msg << "Synopsis must have at least "
                << cell_num_ << "cells"
                << " in the dimension " << syn_dim.getBaseName();
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << err_msg.str();
    }

    // Synopsis chunks info
    chunk_size_ = syn_dim.getChunkInterval();
    chunk_num_ = syn_dim.getCurrEnd() / chunk_size_ + 1;

    // The second synopsis array dimension is the DFT coordinates
    const DimensionDesc &dft_dim = synopsis_desc.getDimensions()[1];
    if (dft_dim.getStartMin() != 0) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Synopsis coordinates must start from 0";
    }
    dft_num_ = dft_dim.getCurrEnd() + 1;

    // Determine low/high attributes
    const Attributes &attrs = synopsis_desc.getAttributes(true);
    if (!SearchArrayDesc::FindAttributeId(attrs, std::string("low"), low_id_) ||
        !SearchArrayDesc::FindAttributeId(attrs, std::string("high"), high_id_)) {
        std::ostringstream err_msg;
        err_msg << "Cannot find low/high attribute in the sample: sample="
                << synopsis_desc.getName();
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << err_msg.str();
    }
}

void Sampler::DFTSynopsis::ParseDFTSynopsisParams(const std::string &params) {
    typedef boost::tokenizer<boost::char_separator<char>> tokenizer_t;
    boost::char_separator<char> sep("xX"); // size_1xsize_2x...xsize_n
    tokenizer_t tokenizer(params, sep);

    // Expecting: AxB, where A is subsequence size and B is MBR size
    tokenizer_t::const_iterator cit = tokenizer.begin();
    if (cit == tokenizer.end()) {
        std::ostringstream err_msg;
        err_msg << "Cannot parse subsequence size from DFT params: " << params;
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << err_msg.str();
    }
    subseq_size_ = boost::lexical_cast<size_t>(cit->c_str());
    ++cit;
    if (cit == tokenizer.end()) {
        std::ostringstream err_msg;
        err_msg << "Cannot parse MBR size from DFT params: " << params;
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << err_msg.str();
    }
    cell_size_ = boost::lexical_cast<size_t>(cit->c_str());
}

void Sampler::DFTSynopsis::SetCacheType(CachingType mode) {
    cache_type_ = mode;
    if (mode == CachingType::EAGER) {
        /*
         * Create the entire cell cache right now to simplify concurrency
         * later. The cells are invalid and will be loaded from the synopsis
         * later in lazy fashion.
         */
        chunks_ = Chunks({true});
        chunks_[0].cells_ = DFTCells(cell_num_);
    } else if (mode == CachingType::LAZY) {
        /*
         * Reserve space for array chunks. All chunks are invalid at the
         * beginning.
         */
        chunks_ = Chunks(chunk_num_);
    } else {
        chunks_ = Chunks(); // make it empty
    }
}

void Sampler::DFTSynopsis::Preload() {
    if (cache_type_ != CachingType::EAGER) {
        LOG4CXX_WARN(logger, "Attempting to preload non-eager synopsis.");
        return;
    }

    // DFT synopsis is linear
    AccessContext ctx;
    for (size_t pos = 0; pos < cell_num_; ++pos) {
    	GetCell(pos, ctx);
    }
    preloaded_ = true;
}

Sampler::DFTCell &Sampler::DFTSynopsis::GetCell(size_t cell_id,
		AccessContext &ctx) {
    if (cache_type_ == CachingType::EAGER) {
        assert(chunks_.size() == 1 && cell_id < chunks_[0].cells_.size());
        DFTCell &cell = chunks_[0].cells_[cell_id];

        // Load cell if needed (cannot check validity here -- concurrency issues)
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
                    FillCellFromArray(cell_id, ctx.iters_, cell);
                    cell.valid_.store(true, std::memory_order_release);
                }
            }
        }
        assert(cell.valid_);
        return cell;
    } else if (cache_type_ == CachingType::LAZY) {
        // Determine position
    	const Coordinate chunk_pos = cell_id - cell_id % chunk_size_;
    	const size_t chunk_linear_pos = chunk_pos / chunk_size_;
        DFTSynopsisChunk &chunk = chunks_[chunk_linear_pos];

        // Fetch (see the double-checked locking comment above)
        bool valid = chunk.valid_.load(std::memory_order_acquire);
        if (!valid) {
            std::lock_guard<std::mutex> lock{mtx_};
            valid = chunk.valid_.load(std::memory_order_relaxed);
            if (!valid) {
                LOG4CXX_DEBUG(logger, "Lazy-loading chunk " << chunk_linear_pos);
                for (size_t i = 0; i < chunk_size_; ++i) {
                	DFTCell &cell = chunk.cells_[i];
                	FillCellFromArray(chunk_pos + i, ctx.iters_, cell);
                	// Relaxed store is fine -- everything is locked
                	cell.valid_.store(true, std::memory_order_relaxed);
                }
                chunk.valid_.store(true, std::memory_order_release);
            }
        }

        // The resulting cell
        assert(chunk.cells_[cell_id - chunk_pos].valid_);
        return chunk.cells_[cell_id - chunk_pos];
    } else {
        FillCellFromArray(cell_id, ctx.iters_, ctx.cell_);
        return ctx.cell_;
    }
}

void Sampler::DFTSynopsis::CheckBounds(Coordinate point) const {
	if (point < synopsis_origin_) {
		point = synopsis_origin_;
	} else if (point > synopsis_end_) {
		point = synopsis_end_;
	}
}

IntervalValue Sampler::DFTSynopsis::SqDist(Coordinate low, Coordinate high,
		const DoubleVector &point) {
	assert(low <= high);
	CheckBounds(low);
	CheckBounds(high);

	// Determine start/end MBRs (synopsis cells)
	const size_t start_cell = (low - synopsis_origin_) / cell_size_;
	const size_t end_cell = (high - synopsis_origin_) / cell_size_;
	AccessContext ctx;
	IntervalValue res;
	res.min_ = std::numeric_limits<double>::max();
	for (size_t cell_id = start_cell; cell_id <= end_cell; ++cell_id) {
		const DFTCell &cell = GetCell(cell_id, ctx);
		// Check for empty MBR
		if (!cell.mbr_.low_.empty()) {
			res.state_ = IntervalValue::NON_NULL;
			const double min_mbr_dist = cell.mbr_.MinSqDist(point);
			if (min_mbr_dist < res.min_) {
				res.min_ = min_mbr_dist;
			}
		}
	}
	return res;
}

void Sampler::DFTSynopsis::FillCellFromArray(size_t pos,
        ArrayIterators &iters, DFTCell &cell) const {

    // Init iterators (first time only)
    if (!iters.low_it_) {
        InitIterators(iters);
    }

    Coordinates syn_pos({Coordinate(pos), 0});
    if (iters.low_it_->setPosition(syn_pos) && !iters.low_it_->isEmpty()) {
    	cell.mbr_.low_.resize(dft_num_);
    	cell.mbr_.high_.resize(dft_num_);
    	for (size_t i = 0; i < dft_num_; ++i) {
    		syn_pos[1] = i;
			// should be a non-empty chunk for sure
			if (!iters.low_it_->setPosition(syn_pos) ||
					!iters.high_it_->setPosition(syn_pos)) {
				std::ostringstream err_msg;
				err_msg << "Cannot get info from synopsis, pos=(";
				std::copy(syn_pos.begin(), syn_pos.end(),
						std::ostream_iterator<Coordinate>(err_msg, ", "));
				err_msg << ")";

				throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
						SCIDB_LE_ILLEGAL_OPERATION) << err_msg.str();
			}
			cell.mbr_.low_[i] = iters.low_it_->getItem().getDouble();
			cell.mbr_.high_[i] = iters.high_it_->getItem().getDouble();
    	}
    }
}

void Sampler::DFTSynopsis::InitIterators(ArrayIterators &iters) const {
    /*
     * One thing to consider here. Creating item iterator results in
     * fetching the first array chunk, which might create a small performance
     * penalty. Another solution is to use array iterators and create
     * chunk iterators when fetching a synopsis cell (see commented code
     * for the FillCell... function).
     */
    iters.low_it_ = synopsis_array_->getItemIterator(low_id_);
    iters.high_it_ = synopsis_array_->getItemIterator(high_id_);
}

size_t Sampler::DFTSynopsis::ComputeMemoryFootprint() const {
    if (cache_type_ == CachingType::EAGER) {
        return cell_num_ * sizeof(DFTCell) + sizeof(DFTSynopsisChunk);
    } else if (cache_type_ == CachingType::LAZY) {
        size_t res = chunk_num_ * sizeof(DFTSynopsisChunk);
        // Count valid chunks' cells.
        for (const auto &c: chunks_) {
            if (c.valid_.load(std::memory_order_relaxed)) {
                res += c.cells_.size() * sizeof(DFTCell);
            }
        }
        return res;
    } else {
        // None
        return 0;
    }
}

Sampler::Synopsis::Synopsis(const ArrayDesc &data_desc,
        const ArrayPtr &array) : synopsis_array_(array) {
    // Array descriptor
    const ArrayDesc &synopsis_desc = array->getArrayDesc();
    const size_t dims_num = synopsis_desc.getDimensions().size();

    // By convenience we store sizes in the name after the last '_'
    const std::string &synopsis_config =
            ArrayDesc::makeUnversionedName(synopsis_desc.getName());
    cell_size_.resize(dims_num);
    ParseChunkSizes(ParseArrayParamsFromName(synopsis_config).back());

    // Compute sample boundaries (in original coordinates)
    synopsis_origin_.resize(dims_num);
    synopsis_end_.resize(dims_num);
    cell_nums_.resize(dims_num);
    chunk_nums_.resize(dims_num);
    chunk_sizes_.resize(dims_num);
    for (size_t i = 0; i < dims_num; i++) {
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
        const uint64_t curr_length = synopsis_end_[i] - synopsis_origin_[i] + 1;
        cell_nums_[i] = (curr_length - 1) / cell_size_[i] + 1;

        // Synopsis correctness checking
        const DimensionDesc &syn_dim = synopsis_desc.getDimensions()[i];
        if (syn_dim.getStartMin() != 0) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Synopsis coordinates must start from 0";
        }
        if (syn_dim.getCurrEnd() + 1 < cell_nums_[i]) { // StartMin() == 0
            std::ostringstream err_msg;
            err_msg << "Synopsis must have at least "
                    << cell_nums_[i] << "cells"
                    << " in the dimension " << syn_dim.getBaseName();
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << err_msg.str();
        }

        // Synopsis chunks info
        chunk_sizes_[i] = syn_dim.getChunkInterval();
        chunk_nums_[i] = syn_dim.getCurrEnd() / chunk_sizes_[i] + 1;
    }

    /*
     * Find min/max/sum/count ids.
     */
    const Attributes &attrs = synopsis_desc.getAttributes(true);
    if (!SearchArrayDesc::FindAttributeId(attrs, std::string("min"), min_id_) ||
        !SearchArrayDesc::FindAttributeId(attrs, std::string("max"), max_id_) ||
        !SearchArrayDesc::FindAttributeId(attrs,
                std::string("count"), count_id_) ||
        !SearchArrayDesc::FindAttributeId(attrs, std::string("sum"), sum_id_)) {
        std::ostringstream err_msg;
        err_msg << "Cannot find min/max attribute in the sample: sample="
                << synopsis_desc.getName();
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << err_msg.str();
    }
}

void Sampler::Synopsis::SetCacheType(CachingType mode) {
    cache_type_ = mode;
    if (mode == CachingType::EAGER) {
        /*
         * Create the entire cell cache right now to simplify concurrency
         * later. The cells are invalid and will be loaded from the synopsis
         * later in lazy fashion.
         */
        chunks_ = Chunks({true});
        chunks_[0].cells_ = Cells(GetTotalCellCount());
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

void Sampler::Synopsis::Preload() {
    if (cache_type_ != CachingType::EAGER) {
        LOG4CXX_WARN(logger, "Attempting to preload non-eager synopsis.");
        return;
    }

    RegionIterator iter{*this, synopsis_origin_, synopsis_end_};
    while (!iter.end()) {
        iter.GetCell();
        ++iter;
    }
    preloaded_ = true;
}

Sampler::Region Sampler::Synopsis::GetSynopsisMBR(const Region &reg) const {
    Region res(reg);
    for (size_t i = 0; i < reg.low_.size(); ++i) {
        res.low_[i] -= (reg.low_[i] - synopsis_origin_[i]) % cell_size_[i];
        res.high_[i] -= (reg.high_[i] - synopsis_origin_[i]) % cell_size_[i];
        res.high_[i] += cell_size_[i] - 1;
    }
    return res;
}

size_t Sampler::Synopsis::GetRegionCost(const Region &reg) const {
    Region mbr(GetSynopsisMBR(reg));

    // MBR is aligned with the synopsis grid, the cost is simple division
    size_t cost = 1;
    for (size_t i = 0; i < mbr.low_.size(); ++i) {
        const size_t len = mbr.high_[i] - mbr.low_[i] + 1;
        assert(len % cell_size_[i] == 0);
        cost *= len / cell_size_[i];
    }
    return cost;
}

Sampler::Cell &Sampler::Synopsis::GetCell(
        const RegionIterator &iter, AccessContext &ctx) {
    if (cache_type_ == CachingType::EAGER) {
        assert(chunks_.size() == 1 &&
                iter.GetCellLinear() < chunks_[0].cells_.size());
        Cell &cell = chunks_[0].cells_[iter.GetCellLinear()];

        // Load cell if needed (cannot check validity here -- concurrency issues)
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
                    FillCellFromArray(iter.GetCurrentSynopsisPosition(),
                            ctx.iters_, cell);
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
                LOG4CXX_DEBUG(logger, "Lazy-loading chunk " << pos.chunk_linear_);
                FillCellsFromChunk(pos, ctx.iters_);
                chunk.valid_.store(true, std::memory_order_release);
            }
        }

        // The cell
        assert(chunk.cells_[pos.cell_chunk_linear_].valid_);
        return chunk.cells_[pos.cell_chunk_linear_];
    } else {
        FillCellFromArray(iter.GetCurrentSynopsisPosition(), ctx.iters_,
                ctx.cell_);
        return ctx.cell_;
    }
}

void Sampler::Synopsis::InitIterators(ArrayIterators &iters) const {
    /*
     * One thing to consider here. Creating item iterator results in
     * fetching the first array chunk, which might create a small performance
     * penalty. Another solution is to use array iterators and create
     * chunk iterators when fetching a synopsis cell (see commented code
     * for the FillCell... function).
     */
    iters.min_it_ = synopsis_array_->getItemIterator(min_id_);
    iters.max_it_ = synopsis_array_->getItemIterator(max_id_);
    iters.sum_it_ = synopsis_array_->getItemIterator(sum_id_);
    iters.count_it_ = synopsis_array_->getItemIterator(count_id_);
}

void Sampler::Synopsis::FillCellsFromChunk(const ChunkCellCoordinates &pos,
        ArrayIterators &iters) const {
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
        FillCellFromArray(chunk_cell_pos.cell_, iters, cell);
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

void Sampler::Synopsis::FillCellFromArray(const Coordinates &pos,
        ArrayIterators &iters, Cell &cell) const {
    // Init iterators (first time only)
    if (!iters.count_it_) {
        InitIterators(iters);
    }

    if (!iters.count_it_->setPosition(pos) || iters.count_it_->isEmpty()) {
        // No chunk in the array -- assume the cell is empty
        cell.count_ = 0;
    } else {
        // should be a non-empty chunk for sure
        if (!iters.min_it_->setPosition(pos) ||
            !iters.max_it_->setPosition(pos) ||
            !iters.sum_it_->setPosition(pos)) {

            std::ostringstream err_msg;
            err_msg << "Cannot get info from synopsis, pos=(";
            std::copy(pos.begin(), pos.end(),
                    std::ostream_iterator<Coordinate>(err_msg, ", "));
            err_msg << ")";

            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                    SCIDB_LE_ILLEGAL_OPERATION) << err_msg.str();
        }

        cell.count_ = iters.count_it_->getItem().getUint64();
        cell.min_ = iters.min_it_->getItem().getDouble();
        cell.max_ = iters.max_it_->getItem().getDouble();
        cell.sum_ = iters.sum_it_->getItem().getDouble();
    }
}

void Sampler::Synopsis::ParseChunkSizes(const std::string &size_param) {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
    boost::char_separator<char> sep("xX"); // size_1xsize_2x...xsize_n
    tokenizer_t tokenizer(size_param, sep);

    shape_cell_size_ = 1;
    int i = 0;
    for (tokenizer_t::const_iterator cit = tokenizer.begin();
            cit != tokenizer.end(); cit++) {
        cell_size_[i] = boost::lexical_cast<Coordinate>(cit->c_str());
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

Sampler::Sampler(const ArrayDesc &data_desc,
        const std::vector<ArrayPtr> &synopsis_arrays,
        const SearchlightConfig &sl_config) :
            data_desc_{data_desc},
            sl_config_(sl_config) {
    // Create catalog of describing synopses
    for (const auto &array: synopsis_arrays) {
        const std::string &array_name =
                ArrayDesc::makeUnversionedName(array->getName());
        const auto &params = ParseArrayParamsFromName(array_name);

        // Attribute name is the second to last component of the name
        const std::string &attr_name = *(params.end() - 2);
        array_synopses_[attr_name].push_back(array);
    }

    // Register default aggregates
    aggrs_["avg"] = AverageSampleAggregate::Create;
    aggrs_["sum"] = SumSampleAggregate::Create;
    aggrs_["min"] = MinMaxSampleAggregate::CreateMin;
    aggrs_["max"] = MinMaxSampleAggregate::CreateMax;

    // Read config params
    cell_thr_ = sl_config_.get("searchlight.sampler.cell_thr", 0.0);
    mbr_thr_ = sl_config_.get("searchlight.sampler.mbr_thr", 0.0);
    cell_limit_ = sl_config_.get("searchlight.sampler.cell_limit", 0);
}

Sampler::~Sampler() {
    std::ostringstream stats_str;
    stats_str << "Stats for synopses follows:\n";
    for (size_t i = 0; i < synopses_.size(); i++) {
        stats_str << "\tStats for attribute " << i << ": \n";
        for (const auto &syn: synopses_[i]) {
            syn->OutputStats(stats_str);
            stats_str << "\n";
        }
    }
    if (logger->isInfoEnabled()) {
        logger->info(stats_str.str());
    }
}

//
// Code for using a pair of array/chunk iterators instead of item iterators.
//
//void Sampler::Synopsis::FillCellFromArray(
//        const Coordinates &pos, Sampler::Cell &cell) {
//    if (!count_it_->setPosition(pos)) {
//        // No chunk in the array -- assume the cell is empty
//        cell = Cell();
//        cell.valid_ = true;
//        return;
//    }
//
//    // should be a non-empty chunk for sure
//    if (!min_it_->setPosition(pos) ||
//        !max_it_->setPosition(pos) ||
//        !sum_it_->setPosition(pos)) {
//
//        std::ostringstream err_msg;
//        err_msg << "Cannot get info from synopsis, pos=(";
//        std::copy(pos.begin(), pos.end(),
//                std::ostream_iterator<std::string>(err_msg, ", "));
//        err_msg << ")";
//
//        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
//                SCIDB_LE_ILLEGAL_OPERATION) << err_msg.str();
//    }
//
//    // now, fetch the chunks
//    const auto &count_chunk_it = count_it_->getChunk().getConstIterator(
//            ConstChunkIterator::IGNORE_EMPTY_CELLS ||
//            ConstChunkIterator::IGNORE_OVERLAPS);
//    if (!count_chunk_it->setPosition(pos) || count_chunk_it->isEmpty()) {
//        // No element in the chunk -- assume empty synopsis cell
//        cell = Cell();
//        cell.valid_ = true;
//        return;
//    }
//
//    // Fetch remaining chunks
//    const auto &min_chunk_it = min_it_->getChunk().getConstIterator(
//            ConstChunkIterator::IGNORE_EMPTY_CELLS ||
//            ConstChunkIterator::IGNORE_OVERLAPS);
//    const auto &max_chunk_it = max_it_->getChunk().getConstIterator(
//            ConstChunkIterator::IGNORE_EMPTY_CELLS ||
//            ConstChunkIterator::IGNORE_OVERLAPS);
//    const auto &sum_chunk_it = sum_it_->getChunk().getConstIterator(
//            ConstChunkIterator::IGNORE_EMPTY_CELLS ||
//            ConstChunkIterator::IGNORE_OVERLAPS);
//    if (!min_chunk_it->setPosition(pos) ||
//        !max_chunk_it->setPosition(pos) ||
//        !sum_chunk_it->setPosition(pos)) {
//
//        std::ostringstream err_msg;
//        err_msg << "Cannot get info from synopsis, pos=(";
//        std::copy(pos.begin(), pos.end(),
//                std::ostream_iterator<std::string>(err_msg, ", "));
//        err_msg << ")";
//
//        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
//                SCIDB_LE_ILLEGAL_OPERATION) << err_msg.str();
//    }
//
//    cell.valid_ = true;
//    cell.count_ = count_chunk_it->getItem().getUint64();
//    cell.min_ = min_chunk_it->getItem().getDouble();
//    cell.max_ = max_chunk_it->getItem().getDouble();
//    cell.sum_ = sum_chunk_it->getItem().getDouble();
//}

StringVector Sampler::ParseArrayParamsFromName(
        const std::string &array_name) {
    StringVector res;

    typedef boost::tokenizer<boost::char_separator<char>> tokenizer_t;
    boost::char_separator<char> sep{"_"}; // parts are separated by '_'
    tokenizer_t tokenizer{array_name, sep};

    for (tokenizer_t::const_iterator cit = tokenizer.begin();
            cit != tokenizer.end(); cit++) {
        res.push_back(*cit);
    }

    if (res.size() < 3) {
        std::ostringstream err_msg;
        err_msg << "Incorrect name for sample array. "
                "Must be name_attr_NxNx...: name=" << array_name;
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << err_msg.str();
    }
    return res;
}

template<class T>
void Sampler::PrepareSynopses(std::vector<T> &synopses,
		CachingType cache_type, bool preload, size_t mem_limit,
		const std::string &attr_name, AttributeID attr_search_id) {
	// First, sort synopses by their cell sizes
    std::stable_sort(synopses.begin(), synopses.end(),
            [](const T &s1, const T &s2) {
        return s1->GetCellSize() > s2->GetCellSize();
    });

    // Then set cache type and preload
	for (const auto &syn: synopses) {
		const size_t syn_mem_mb = syn->MemorySize() / 1024 / 1024; // in MB
		if (syn_mem_mb <= mem_limit) {
			// Cache the synopsis
			syn->SetCacheType(cache_type);
			// .. and preload it if needed
			if (preload) {
				LOG4CXX_INFO(logger, "Preloading synopsis: "<< syn->GetName());
				syn->Preload();
			}
			mem_limit -= syn_mem_mb;
		} else {
			break;
		}
	}

    // Warn about performance problems...
    if (!synopses.front()->IsCached()) {
        LOG4CXX_WARN(logger, "No synopses are cached for attribute "
                << attr_name << "(" << attr_search_id << ")");
    }

    // Debug printing
    if (logger->isDebugEnabled()) {
        std::ostringstream msg;
        msg << "Synopses loaded for attribute " << attr_name
                << '(' << attr_search_id << "): ";
        for (const auto &syn: synopses) {
            msg << syn->GetName() << "(cached: " << int(syn->GetCachingType())
                    << "), ";
        }
        logger->debug(msg.str());
    }
}

void Sampler::LoadSampleForAttribute(const std::string &attr_name,
        AttributeID attr_search_id) {
    // Add a new vector of synopses if needed
    if (attr_search_id + 1 > synopses_.size()) {
        synopses_.resize(attr_search_id + 1);
        dft_synopses_.resize(attr_search_id + 1);
    }

    // If we have loaded something for the attribute, ignore
    auto &loaded_synopses = synopses_[attr_search_id];
    if (!loaded_synopses.empty()) {
        return;
    }

    // See what we have for the attribute
    const auto &synops = array_synopses_[attr_name];
    auto &loaded_dft_synopses = dft_synopses_[attr_search_id];
    for (const auto &syn: synops) {
    	// Determine if we're loading a DFT synopsis: they start as "dft_"
        const std::string &syn_name =
                ArrayDesc::makeUnversionedName(syn->getName());
    	if (ParseArrayParamsFromName(syn_name)[0] == "dft") {
    		loaded_dft_synopses.emplace_back(new DFTSynopsis{data_desc_, syn});
    	} else {
    		loaded_synopses.emplace_back(new Synopsis{data_desc_, syn});
    	}
    }

    // Cache type
    const std::string cache_type_string =
            sl_config_.get("searchlight.sampler.cache", "eager");
    CachingType cache_type = CachingType::NONE;
    if (cache_type_string == "eager") {
        cache_type = CachingType::EAGER;
    } else if (cache_type_string == "lazy") {
        cache_type = CachingType::LAZY;
    } else {
        LOG4CXX_WARN(logger, "Unknown cache type: " << cache_type_string);
    }

	// Are we preloading synopses?
    const bool preload_syns =
			sl_config_.get("searchlight.sampler.preload", 0);

    // Set cache type and preload aggregate synopses
    if (!loaded_synopses.empty()) {
		const size_t memory_limit_mb =
				sl_config_.get("searchlight.sampler.memory_per_attr", 1024);
    	PrepareSynopses(loaded_synopses, cache_type, preload_syns,
    			memory_limit_mb, attr_name, attr_search_id);
    }

    // Set cache type and preload for DFT synopses
    if (!loaded_dft_synopses.empty()) {
		const size_t memory_limit_mb =
				sl_config_.get("searchlight.sampler.memory_per_attr_dft", 1024);
    	PrepareSynopses(loaded_dft_synopses, cache_type, preload_syns,
    			memory_limit_mb, attr_name, attr_search_id);
    }
}

size_t Sampler::Synopsis::ComputeMemoryFootprint() const {
    if (cache_type_ == CachingType::EAGER) {
        return GetTotalCellCount() * sizeof(Cell) + sizeof(SynopsisChunk);
    } else if (cache_type_ == CachingType::LAZY) {
        size_t res = GetTotalChunksCount() * sizeof(SynopsisChunk);
        // Count valid chunks' cells.
        for (const auto &c: chunks_) {
            if (c.valid_.load(std::memory_order_relaxed)) {
                res += c.cells_.size() * sizeof(Cell);
            }
        }
        return res;
    } else {
        // None
        return 0;
    }
}

void Sampler::Synopsis::CheckAndCorrectBounds(
        Coordinates &low, Coordinates &high) const {
    if (low.size() != high.size()) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Specified region has inconsistent dimensions!";
    }

    for (size_t i = 0; i < low.size(); i++) {
        if (low[i] > high[i]) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Specified region has low > high coordinates!";
        }
        if (low[i] < synopsis_origin_[i]) {
            low[i] = synopsis_origin_[i];
        }
        if (high[i] > synopsis_end_[i]) {
            high[i] = synopsis_end_[i];
        }
    }
}

bool Sampler::Synopsis::CheckBounds(const Coordinates &point) const {
    if (point.size() != synopsis_origin_.size()) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Specified point has inconsistent dimensions!";
    }

    for (size_t i = 0; i < synopsis_origin_.size(); i++) {
        if (point[i] < synopsis_origin_[i] || point[i] > synopsis_end_[i]) {
            return false;
        }
    }

    return true;
}

bool Sampler::Synopsis::RegionValidForSample(const Coordinates &low,
        const Coordinates &high) const {
    // copy to check correctness and align
    Coordinates inlow(low), inhigh(high);
    CheckAndCorrectBounds(inlow, inhigh);

    /*
     * The synopsis can handle the region if:
     *   1) Its corner (low) is aligned with a synopsis cell.
     *   2) Each side length is divisible by the cell size.
     *
     * That means the region contains only full synopsis cells, no partial
     * intersections.
     */
    for (size_t i = 0; i < inlow.size(); i++) {
        const Coordinate cs = cell_size_[i];
        const Coordinate len = inhigh[i] - inlow[i] + 1;
        const Coordinate low_off = inlow[i] - synopsis_origin_[i];
        if (low_off % cs != 0 || len % cs != 0) {
            return false;
        }
    }
    return true;
}

void Sampler::Synopsis::ComputeAggregate(const Coordinates &low,
        const Coordinates &high, const SampleAggregatePtrVector &aggs) {
    // copy to check correctness and align
    Coordinates inlow(low), inhigh(high);
    CheckAndCorrectBounds(inlow, inhigh);

    // go through the chunks
    RegionIterator iter(*this, inlow, inhigh);
    uint64_t cells_accessed = 0;
    while (!iter.end()) {
        const Cell &cell = iter.GetCell();
        cells_accessed++;
        const uint64_t part_size = iter.GetPartSize();
        for (size_t i = 0; i < aggs.size(); i++) {
            if (cell.count_ > 0) {
                aggs[i].get()->AccumulateChunk(shape_cell_size_,
                        part_size, cell);
            }
        }
        ++iter;
    }
    cells_accessed_.fetch_add(cells_accessed, std::memory_order_relaxed);
}

void Sampler::Synopsis::ComputeAggregatesWithThr(
        const SampleAggregatePtrVector &aggs,
        const std::vector<Region> &in_regions,
        std::vector<Region> &left_regions, double cell_thr) {
    // We will reuse the iterator, but have to init it first with some values
    RegionIterator iter{*this, synopsis_origin_, synopsis_origin_};
    uint64_t cells_accessed = 0;
    for (const auto &region: in_regions) {
        iter.Reset(region.low_, region.high_);
        while (!iter.end()) {
            const uint64_t part_size = iter.GetPartSize();
            if (double(part_size) / shape_cell_size_ >= cell_thr) {
                const Cell &cell = iter.GetCell();
                cells_accessed++;
                for (size_t i = 0; i < aggs.size(); i++) {
                    if (cell.count_ > 0) {
                        aggs[i].get()->AccumulateChunk(shape_cell_size_,
                                part_size, cell);
                    }
                }
            } else {
                Region part;
                iter.GetCoveredCell(part.low_, part.high_);
                left_regions.push_back(std::move(part));
            }
            ++iter;
        }
    }
    cells_accessed_.fetch_add(cells_accessed, std::memory_order_relaxed);
}

IntervalValue Sampler::Synopsis::GetElement(const Coordinates &point) {
    IntervalValue res;
    if (!CheckBounds(point)) {
        return res; // out of bounds -- NULL
    }

    RegionIterator iter(*this, point, point);
    const Cell &cell = iter.GetCell();
    cells_accessed_.fetch_add(1, std::memory_order_relaxed);

    if (cell.count_ == 0) {
        return res; // in an empty chunk -- definitely NULL
    } else if (cell.count_ == shape_cell_size_) {
        res.state_ = IntervalValue::NON_NULL; // full chunk -- no empty elems
    } else {
        res.state_ = IntervalValue::MAY_NULL;
    }

    res.max_ = cell.max_;
    res.min_ = cell.min_;
    res.val_ = cell.sum_ / cell.count_; // uniformity assumption

    return res;
}

IntervalValueVector Sampler::ComputeAggregate(const Coordinates &low,
        const Coordinates &high, AttributeID s_attr,
        const StringVector &aggr_names, bool exact) const {
    assert(s_attr < synopses_.size());

    // resolve aggregates
    SampleAggregatePtrVector aggs(aggr_names.size());
    for (size_t i = 0; i < aggr_names.size(); i++) {
        AggrMap::const_iterator it = aggrs_.find(aggr_names[i]);
        if (it == aggrs_.end()) {
            std::ostringstream err_msg;
            err_msg << "Sample aggregate not found: aggr=" << aggr_names[i];
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << err_msg.str();
        }
        // the aggregate is created via the registered factory
        aggs[i].reset(it->second());
    }

    // Choose synopsis
    Synopsis *syn = nullptr;
    if (exact) {
        for (const auto &s: synopses_[s_attr]) {
            if (s->RegionValidForSample(low, high)) {
                syn = s.get();
                break;
            }
        }
        /*
         * If we cannot find a suitable synopsis, that means the exact
         * computation cannot be performed by the sampler.
         */
        if (!syn) {
            return IntervalValueVector();
        } else {
            syn->ComputeAggregate(low, high, aggs);
        }
    } else {
        // Non-exact computation: use the threshold algorithm
        assert(!synopses_[s_attr].empty());
        auto &syns = synopses_[s_attr];

        const Region orig_region{low, high};
        std::vector<Region> regions{orig_region};
        syns[0]->CheckAndCorrectBounds(regions[0].low_, regions[0].high_);
        std::vector<Region> left_regions;

        for (size_t i = 0; i < syns.size(); i++) {
            // Check if we have a cached higher-resolution synopsis
            bool have_better_synopses = (i != syns.size() - 1) &&
                    syns[i + 1]->IsCached();

            // Check the MBR threshold
            if (mbr_thr_ > 0.0) {
                const Region syn_mbr(syns[i]->GetSynopsisMBR(orig_region));
                if (orig_region.AreaRatio(syn_mbr) >= mbr_thr_) {
                    have_better_synopses = false;
                }
            }

            syns[i]->ComputeAggregatesWithThr(aggs, regions, left_regions,
                    have_better_synopses ? cell_thr_ : 0.0);
            if (!left_regions.empty()) {
                regions = left_regions;
                left_regions.clear();
            } else {
                break;
            }
        }
        assert(left_regions.empty());
    }

    // finalize the result
    IntervalValueVector res(aggs.size());
    for (size_t i = 0; i < aggs.size(); i++) {
        aggs[i].get()->Finalize(res[i]);
    }

    return res;
}

IntervalValue Sampler::GetElement(const Coordinates &point,
        AttributeID attr) const {
    // We always estimate element via the primary sample
    assert(!synopses_[attr].empty());
    const auto &syn = synopses_[attr][0];
    return syn->GetElement(point);
}

IntervalValue Sampler::SqDist(AttributeID attr, Coordinate low, Coordinate high,
		const DoubleVector &point) const {
    // FIXME: For now we compute distance via the primary sample
	assert(!dft_synopses_[attr].empty());
	const auto &syn = dft_synopses_[attr][0];
	return syn->SqDist(low, high, point);
}
} /* namespace searchlight */

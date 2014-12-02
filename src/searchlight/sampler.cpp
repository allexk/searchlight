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
    for (size_t i = 0; i < dims_num; i++) {
        // Metadata is filled from the data descriptor
        const DimensionDesc &data_dim = data_desc.getDimensions()[i];
        if (data_dim.getStartMin() == scidb::MIN_COORDINATE ||
                data_dim.getEndMax() == scidb::MAX_COORDINATE) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Unbounded arrays are not supported by the sampler!";
        }
        synopsis_origin_[i] = data_dim.getStartMin();
        synopsis_end_[i] = data_dim.getEndMax();
        cell_nums_[i] = (data_dim.getLength() - 1) / cell_size_[i] + 1;

        // Synopsis correctness checking
        const DimensionDesc &dim = synopsis_desc.getDimensions()[i];
        if (dim.getStartMin() != 0) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Synopsis coordinates must start from 0";
        }
        if (dim.getLength() != cell_nums_[i]) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Synopsis must have " << cell_nums_[i] << "cells"
                    << " in the dimension " << dim.getBaseName();
        }
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

void Sampler::Synopsis::SetCacheMode(bool mode) {
    cache_cells_ = mode;
    if (mode) {
        /*
         * Create the entire cell cache right now to simplify concurrency
         * later. The cells are invalid and will be loaded from the synopsis
         * later in lazy fashion.
         */
        cells_.clear();
        cells_.resize(GetTotalCellCount());
    } else {
        cells_.clear();
        cells_.shrink_to_fit();
    }
}

void Sampler::Synopsis::Preload() {
    if (!cache_cells_) {
        LOG4CXX_WARN(logger, "Attempting to preload non-cached synopsis.");
        return;
    }

    RegionIterator iter{*this, synopsis_origin_, synopsis_end_};
    while (!iter.end()) {
        iter.GetCell();
        ++iter;
    }
    preloaded_ = true;
}

const Sampler::Cell &Sampler::Synopsis::CachedAccessor::GetCell(
        const RegionIterator &iter) {
    return syn_.CacheCell(iter, iters_);
}

const Sampler::Cell &Sampler::Synopsis::NonCachedAccessor::GetCell(
        const RegionIterator &iter) {
    syn_.FillCellFromArray(iter.GetCurrentSynopsisPosition(), iters_, cell_);
    return cell_;
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

void Sampler::Synopsis::LoadCellFromArray(const Coordinates &pos,
        ArrayIterators &iters, Cell &cell) const {
    /*
     * Need a lock here. For efficiency reasons this has been
     * implemented as a double-checked locking with fences and atomics.
     */
    bool valid = cell.valid_.load(std::memory_order_acquire);
    if (!valid) {
        std::lock_guard<std::mutex> lock{mtx_};
        valid = cell.valid_.load(std::memory_order_relaxed);
        if (!valid) {
            FillCellFromArray(pos, iters, cell);
            cell.valid_.store(true, std::memory_order_release);
        }
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

void Sampler::LoadSampleForAttribute(const std::string &attr_name,
        AttributeID attr_search_id) {
    // Add a new vector of synopses if needed
    if (attr_search_id + 1 > synopses_.size()) {
        synopses_.resize(attr_search_id + 1);
    }

    // If we have loaded something for the attribute, ignore
    auto &loaded_synopses = synopses_[attr_search_id];
    if (!loaded_synopses.empty()) {
        return;
    }

    // See what we have for the attribute
    const auto &synops = array_synopses_[attr_name];
    for (const auto &syn: synops) {
        loaded_synopses.emplace_back(new Synopsis{data_desc_, syn});
    }

    if (!loaded_synopses.empty()) {
        /*
         *  Sort the vector of synopses by the cell size, preserving user's
         *  ordering for synopses of the same cell size.
         */
        std::stable_sort(loaded_synopses.begin(), loaded_synopses.end(),
                [](const SynopsisPtr &s1, const SynopsisPtr &s2) {
            return s1->GetCellSize() > s2->GetCellSize();
        });

        // Enable caching for some of the synopses
        size_t memory_limit_mb =
                sl_config_.get("searchlight.sampler.memory_per_attr", 1024);
        const bool preload_syns =
                sl_config_.get("searchlight.sampler.preload", 1);
        for (const auto &syn: loaded_synopses) {
            const size_t syn_mem_mb = syn->MemorySize() / 1024 / 1024; // in MB
            if (syn_mem_mb <= memory_limit_mb) {
                // Cache the synopsis
                syn->SetCacheMode(true);
                // .. and preload it if needed
                if (preload_syns) {
                    LOG4CXX_INFO(logger, "Preloading synopsis: "
                            << syn->GetName());
                    syn->Preload();
                }
                memory_limit_mb -= syn_mem_mb;
            } else {
                break;
            }
        }

        // Warn about performance problems...
        if (!loaded_synopses.front()->IsCached()) {
            LOG4CXX_WARN(logger, "No synopses are cached for attribute "
                    << attr_name << "(" << attr_search_id << ")");
        }

        // Debug printing
        if (logger->isDebugEnabled()) {
            std::ostringstream msg;
            msg << "Synopses loaded for attribute " << attr_name
                    << '(' << attr_search_id << "): ";
            for (const auto &syn: loaded_synopses) {
                msg << syn->GetName() << "(cached: " << syn->IsCached()
                        << "), ";
            }
            logger->debug(msg.str());
        }
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

const Sampler::Cell &Sampler::Synopsis::CacheCell(
        const RegionIterator &iter, ArrayIterators &iters) {
    assert(iter.GetCellLinear() < cells_.size());
    Cell &cell = cells_[iter.GetCellLinear()];

    // Load cell if needed (cannot check validity here -- concurrency issues)
    if (!preloaded_) {
        LoadCellFromArray(iter.GetCurrentSynopsisPosition(), iters, cell);
    }
    assert(cell.valid_);

    return cell;
}

IntervalValueVector Sampler::Synopsis::ComputeAggregate(const Coordinates &low,
        const Coordinates &high, const SampleAggregatePtrVector &aggs) {
    // copy to check correctness and align
    Coordinates inlow(low), inhigh(high);
    CheckAndCorrectBounds(inlow, inhigh);

    // go through the chunks
    RegionIterator iter(*this, inlow, inhigh);
    while (!iter.end()) {
        const Cell &cell = iter.GetCell();
        const uint64_t part_size = iter.GetPartSize();
        for (size_t i = 0; i < aggs.size(); i++) {
            if (cell.count_ > 0) {
                aggs[i].get()->AccumulateChunk(shape_cell_size_,
                        part_size, cell);
            }
        }
        ++iter;
    }

    // finalize the result
    IntervalValueVector res(aggs.size());
    for (size_t i = 0; i < aggs.size(); i++) {
        aggs[i].get()->Finalize(res[i]);
    }

    return res;
}

IntervalValue Sampler::Synopsis::GetElement(const Coordinates &point) {
    IntervalValue res;
    if (!CheckBounds(point)) {
        return res; // out of bounds -- NULL
    }

    RegionIterator iter(*this, point, point);
    const Cell &cell = iter.GetCell();

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
        }
    } else {
        // For non-exact computation alway use the primary synopsis
        assert(!synopses_[s_attr].empty());
        syn = synopses_[s_attr][0].get();
    }

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

    return syn->ComputeAggregate(low, high, aggs);
}

IntervalValue Sampler::GetElement(const Coordinates &point,
        AttributeID attr) const {
    // We always estimate element via the primary sample
    assert(!synopses_[attr].empty());
    const auto &syn = synopses_[attr][0];
    return syn->GetElement(point);
}
} /* namespace searchlight */

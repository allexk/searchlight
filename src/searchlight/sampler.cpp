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

class AverageSampleAggregate : public SampleAggregate {
public:
    static SampleAggregate *Create() {
        return new AverageSampleAggregate;
    }

    virtual void AccumulateChunk(uint64_t chunk_size, uint64_t part_size,
            const Sampler::Chunk &chunk) {
        // chunk info
        const double sum = chunk.sum_;
        const uint64_t count = chunk.count_;
        const double min = chunk.min_;
        const double max = chunk.max_;

        // is it a full chunk? guaranteed to be non-empty
        if (part_size == chunk_size) {
            not_null_ = true;
            sum_ += sum;
            approx_sum_ += sum;
            count_ += count;
            approx_count_ += count;
            return;
        }

        // a part of a chunk; guaranteed non-empty
        const size_t chunk_num = chunk_min_count_.size();
        const uint64_t empty_count = chunk_size - count;

        // compute counts
        const uint64_t min_count = part_size >= empty_count ?
                part_size - empty_count : 0;
        const uint64_t max_count = part_size >= count ? count : part_size;
        if (min_count > 0) {
            not_null_ = true;
        }
        chunk_min_count_.push_back(min_count);
        chunk_max_count_.push_back(max_count);

        /*
         * For the approximate computation we make the uniformity assumption,
         * including the uniform spread of empty elements. In practice,
         * it is probably less probable -- empty elements tend to stick
         * together. But without any additional knowledge we cannot do better.
         */
        const double part_ratio = double(part_size) / chunk_size;
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
            chunk_efforts_upper_.push_back(Effort(max, max_num, chunk_num));
            count_left -= max_num;
            const double deltau = sum - max_num * max -
                    (count_left - 1) * min;
            chunk_efforts_upper_.push_back(Effort(deltau, 1, chunk_num));
            count_left--;
            if (count_left > 0) {
                chunk_efforts_upper_.push_back(Effort(min, count_left,
                        chunk_num));
            }

            // lower boundary efforts
            count_left = count;
            const uint64_t min_num = floor((count * max - sum) / (max - min));
            chunk_efforts_lower_.push_back(Effort(min, min_num, chunk_num));
            count_left -= min_num;
            const double deltal = sum - min_num * min -
                    (count_left - 1) * max;
            chunk_efforts_lower_.push_back(Effort(deltal, 1, chunk_num));
            count_left--;
            if (count_left > 0) {
                chunk_efforts_lower_.push_back(Effort(max, count_left,
                        chunk_num));
            }
        } else {
            chunk_efforts_upper_.push_back(Effort(max, count, chunk_num));
            chunk_efforts_lower_.push_back(Effort(min, count, chunk_num));
        }
    }

    virtual void Finalize(IntervalValue &res) {
        if (not_null_) {
            res.state_ = IntervalValue::NON_NULL;
        } else {
            if (!chunk_min_count_.empty()) {
                // got possible non-empty parts
                res.state_ = IntervalValue::MAY_NULL;
            } else {
                res.state_ = IntervalValue::NUL;
                return;
            }
        }

        // approximate value
        res.val_ = approx_sum_ / approx_count_;

        // we have to make copies of counts to reuse them for the lower bound
        std::vector<uint64_t> chunk_min_count_copy(chunk_min_count_);
        std::vector<uint64_t> chunk_max_count_copy(chunk_max_count_);

        // sort the efforts by the value
        std::sort(chunk_efforts_upper_.begin(), chunk_efforts_upper_.end(),
                Effort::GreaterOrder()); // descending
        std::sort(chunk_efforts_lower_.begin(), chunk_efforts_lower_.end());

        /*
         * Compute the upper bound: raise the average while you can, and then
         * lower as less as possible.
         */
        double current_sum = sum_;
        double current_count = count_;
        for (size_t i = 0; i < chunk_efforts_upper_.size(); i++) {
            const Effort &eff = chunk_efforts_upper_[i];
            uint64_t &min_count_chunk = chunk_min_count_[eff.chunk_];
            uint64_t &max_count_chunk = chunk_max_count_[eff.chunk_];
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
        for (size_t i = 0; i < chunk_efforts_lower_.size(); i++) {
            const Effort &eff = chunk_efforts_lower_[i];
            uint64_t &min_count_chunk = chunk_min_count_copy[eff.chunk_];
            uint64_t &max_count_chunk = chunk_max_count_copy[eff.chunk_];
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
     *  An effort describes the chunk's best effort to bring up/down the avg.
     *  For example, to increase the avg. maximally we would add the value_
     *  the number_ times to the set of values composing the upper bound.
     *  Chunk_ references the chunk number for this effort.
     */
    struct Effort {
        double value_;
        uint64_t number_;
        size_t chunk_;

        Effort(double value, uint64_t number, size_t chunk) :
            value_(value), number_(number_), chunk_(chunk) {}

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
    std::vector<Effort> chunk_efforts_upper_, chunk_efforts_lower_;

    /*
     * The minimum and maximum number of elements for each chunk.
     */
    std::vector<uint64_t> chunk_min_count_, chunk_max_count_;

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

    virtual void AccumulateChunk(uint64_t chunk_size,
            uint64_t part_size, const Sampler::Chunk &chunk) {
        // chunk info
        const double sum = chunk.sum_;
        const uint64_t count = chunk.count_;
        const double min = chunk.min_;
        const double max = chunk.max_;

        // it cannot be null, since chunks are non-empty here
        null_ = false;

        // is it a full chunk? guaranteed to be non-empty
        if (part_size == chunk_size) {
            not_null_ = true;
            min_sum_ += sum;
            max_sum_ += sum;
            approx_sum_ += sum;
            return;
        }

        // compute counts
        const uint64_t empty_count = chunk_size - count;
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
        const double part_ratio = double(part_size) / chunk_size;
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
        res.val_ = approx_sum_;
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

    virtual void AccumulateChunk(uint64_t chunk_size,
            uint64_t part_size, const Sampler::Chunk &chunk) {
        // chunk info
        const double min = chunk.min_;
        const double max = chunk.max_;

        // it cannot be null, since chunks are non-empty here
        null_ = false;

        // is it a full chunk? guaranteed to be non-empty
        if (part_size == chunk_size) {
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
        const uint64_t empty_count = chunk_size - chunk.count_;
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
        min_(std::numeric_limits<double>::max()),
        max_(std::numeric_limits<double>::min()),
        approx_(0),
        is_min_(min), not_null_(false), null_(true) {}

    // Minimum/maximum and approximate values
    double min_, max_, approx_;

    // is it the min aggregate?
    bool is_min_;

    // is it definitely not null? or is it definitely null?
    bool not_null_, null_;
};

Sampler::Sampler(const ArrayPtr &array, const ArrayDesc &data_desc) :
        sample_array_(array),
        chunk_sizes_(data_desc.getDimensions().size()),
        sample_origin_(data_desc.getDimensions().size()),
        sample_end_(data_desc.getDimensions().size()),
        chunk_nums_(data_desc.getDimensions().size()) {
    const ArrayDesc &sample_desc = array->getArrayDesc();

    // by convenience we store sizes in the comment :)
    const std::string &sample_config =
            ArrayDesc::makeUnversionedName(sample_desc.getName());
    const size_t undersc_pos = sample_config.find_last_of('_');
    if (undersc_pos == std::string::npos || sample_config.size() < 3 ||
            undersc_pos > sample_config.size() - 2) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Incorrect sample array name (must have an underscore "
                        "followed by ?x? size parameters: name="
                << sample_desc.getName();
    }
    SetChunkSizes(sample_config.substr(undersc_pos + 1));

    // The start of the sample corresponds to the start of the data array
    for (size_t i = 0; i < data_desc.getDimensions().size(); i++) {
        const DimensionDesc &dim = data_desc.getDimensions()[i];
        sample_origin_[i] = dim.getCurrStart();
        sample_end_[i] = dim.getCurrEnd();
        chunk_nums_[i] = (sample_end_[i] - sample_origin_[i]) /
                chunk_sizes_[i] + 1;
    }

    // number of elements in a chunk
    shape_chunk_size_ = ComputeShapeChunkSize();

    // find min/max ids (not necessary, unless the array has the empty bitmap)
    const Attributes &attrs = sample_desc.getAttributes(true);
    if (!SearchArrayDesc::FindAttributeId(attrs, std::string("min"), min_id_) ||
        !SearchArrayDesc::FindAttributeId(attrs, std::string("max"), max_id_) ||
        !SearchArrayDesc::FindAttributeId(attrs,
                std::string("count"), count_id_) ||
        !SearchArrayDesc::FindAttributeId(attrs, std::string("sum"), sum_id_)) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Cannot find min/max attribute in the sample: sample="
                << sample_desc.getName();
    }

    if (sample_desc.getDimensions()[0].getCurrStart() != 0) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Chunk coordinate should start from 0: sample="
                << sample_desc.getName();
    }
    chunks_num_ = sample_desc.getDimensions()[0].getCurrEnd();

    // Register default aggregates
    aggrs_["avg"] = AverageSampleAggregate::Create;
    aggrs_["sum"] = SumSampleAggregate::Create;
    aggrs_["min"] = MinMaxSampleAggregate::CreateMin;
    aggrs_["max"] = MinMaxSampleAggregate::CreateMax;
}

void Sampler::LoadSampleForAttribute(AttributeID attr_orig_id,
        AttributeID attr_search_id) {

    boost::shared_ptr<ConstItemIterator> min_iterator =
            sample_array_->getItemIterator(min_id_);
    boost::shared_ptr<ConstItemIterator> max_iterator =
            sample_array_->getItemIterator(max_id_);
    boost::shared_ptr<ConstItemIterator> count_iterator =
            sample_array_->getItemIterator(count_id_);
    boost::shared_ptr<ConstItemIterator> sum_iterator =
            sample_array_->getItemIterator(sum_id_);

    // Sample: first dimension -- region, second -- the original attribute
    sample_chunks_.push_back(ChunkVector());
    if (sample_chunks_.size() != attr_search_id) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Sampler and descriptor inconsistency: sample aid="
                << sample_chunks_.size() << ", desc id=" << attr_search_id;
    }
    ChunkVector &chunks = sample_chunks_.back();
    chunks.reserve(chunks_num_);
    Coordinates pos(2);
    pos[1] = attr_orig_id;
    for (pos[0] = 0; pos[0] < chunks_num_; pos[0]++) {
        if (!min_iterator->setPosition(pos) ||
                !max_iterator->setPosition(pos) ||
                !count_iterator->setPosition(pos) ||
                !sum_iterator->setPosition(pos)) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Cannot get info from sample, chunk=" << pos[0] <<
                    "attr=" << pos[1];
        }
        const double minv = min_iterator->getItem().getDouble();
        const double maxv = max_iterator->getItem().getDouble();
        const double sumv = sum_iterator->getItem().getDouble();
        const uint64_t count = count_iterator->getItem().getUint64();
        chunks.push_back(Chunk(minv, maxv, sumv, count));
    }
}

void Sampler::SetChunkSizes(const std::string &size_param) {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
    boost::char_separator<char> sep("xX"); // size_1xsize_2x...xsize_n
    tokenizer_t tokenizer(size_param, sep);

    int i = 0;
    for (tokenizer_t::const_iterator cit = tokenizer.begin();
            cit != tokenizer.end(); cit++) {
        chunk_sizes_[i++] = boost::lexical_cast<Coordinate>(cit->c_str());
    }

    if (i != chunk_sizes_.size()) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Could not retrieve all chunk sizes: conf=" <<
                size_param << ", needed sizes=" << chunk_sizes_.size();
    }
}

void Sampler::CheckAndCorrectBounds(Coordinates &low, Coordinates &high) const {
    if (low.size() != high.size()) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Specified region has inconsistent dimensions!";
    }

    for (size_t i = 0; i < low.size(); i++) {
        if (low[i] > high[i]) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Specified region has low > high coordinates!";
        }
        if (low[i] < sample_origin_[i]) {
            low[i] = sample_origin_[i];
        }
        if (high[i] > sample_end_[i]) {
            high[i] = sample_end_[i];
        }
    }
}

bool Sampler::CheckBounds(const Coordinates &point) const {
    if (point.size() != sample_origin_.size()) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Specified point has inconsistent dimensions!";
    }

    for (size_t i = 0; i < sample_origin_.size(); i++) {
        if (point[i] < sample_origin_[i] || point[i] > sample_end_[i]) {
            return false;
        }
    }

    return true;
}

uint64_t Sampler::ComputeShapeChunkSize() const {
    uint64_t res = 1;
    for (size_t i = 0; i < chunk_sizes_.size(); i++) {
        res *= chunk_sizes_[i];
    }
    return res;
}

IntervalValueVector Sampler::ComputeAggregate(const Coordinates &low,
        const Coordinates &high, AttributeID attr,
        const StringVector &aggr_names) const {
    // copy to check correctness and align
    Coordinates inlow(low), inhigh(high);
    CheckAndCorrectBounds(inlow, inhigh);

    // resolve aggregates
    SampleAggregatePtrVector aggs(aggr_names.size());
    for (size_t i = 0; i < aggr_names.size(); i++) {
        AggrMap::const_iterator it = aggrs_.find(aggr_names[i]);
        if (it == aggrs_.end()) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Sample aggregate not found: aggr=" << aggr_names[i];
        }
        // the aggregate is created via the registered factory
        aggs[i].reset(it->second());
    }

    // go through the chunks
    RegionIterator iter(*this, inlow, inhigh);
    const ChunkVector &chunks = sample_chunks_[attr];
    while (!iter.end()) {
        const Chunk &chunk = chunks[iter.GetChunk()];
        const uint64_t part_size = iter.GetPartSize();
        for (size_t i = 0; i < aggs.size(); i++) {
            if (chunk.count_ > 0) {
                aggs[i].get()->AccumulateChunk(shape_chunk_size_,
                        part_size, chunk);
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

IntervalValue Sampler::GetElement(const Coordinates &point,
        AttributeID attr) const {
    IntervalValue res;
    if (!CheckBounds(point)) {
        return res; // out of bounds -- NULL
    }

    RegionIterator iter(*this, point, point);
    const Chunk &chunk = sample_chunks_[attr][iter.GetChunk()];

    if (chunk.empty()) {
        return res; // in an empty chunk -- definitely NULL
    } else if (chunk.count_ == shape_chunk_size_) {
        res.state_ = IntervalValue::NON_NULL; // full chunk -- no empty elems
    } else {
        res.state_ = IntervalValue::MAY_NULL;
    }

    res.max_ = chunk.max_;
    res.min_ = chunk.min_;
    res.val_ = chunk.sum_ / chunk.count_; // uniformity assumption

    return res;
}
} /* namespace searchlight */

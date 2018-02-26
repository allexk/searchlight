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

#include <fftw3.h>

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
            const AggCell &cell) override{
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
            uint64_t part_size, const AggCell &cell) override {
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
            uint64_t part_size, const AggCell &cell) override {
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

SeqSynopsis::SeqSynopsis(const ArrayDesc &data_desc,
                         const ArrayPtr &array) :
        Base{data_desc, array,
             array->getArrayDesc().getDimensions().size() - 1} {
    // Array descriptor
    const ArrayDesc &synopsis_desc = array->getArrayDesc();
    // Read the seq-specific config (the Base has already read some)
    const auto synopsis_params = TokenizeString(
            ArrayDesc::makeUnversionedName(synopsis_desc.getName()), "_");
    if (synopsis_params.size() < 3) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION) <<
                "DFT synopsis array must have name:"
                "<type>_<subseq_size>_<name>_<attr>_<cell_params>";
    }
    // For sequence synopsis, type is the first component
    const std::string type_param = synopsis_params[0];
    if (type_param == "dft") {
        type_ = Type::DFT;
    } else if (type_param == "paa") {
        type_ = Type::PAA;
    } else {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION) <<
                "Unknown sequence synopsis type";
    }
    // The second to last contains the subsequence length
    try {
        subseq_size_ = std::stoi(synopsis_params[1]);
    } catch (const std::invalid_argument &) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION) <<
                "The second to last argument must be subsequence size";
    }
    // The remaining parameters are derived from the descriptor
    mbr_size_ = cell_size_.back();
    features_num_ = synopsis_desc.getDimensions().back().getLength();
    // Correct cell_nums_ for the MBR dimension (-2 since the last is "coord")
    const size_t mbr_dim_index = synopsis_desc.getDimensions().size() - 2;
    //const DimensionDesc &mbr_dim = synopsis_desc.getDimensions()[mbr_dim_index];
    const size_t mbr_dim_len = (synopsis_end_[mbr_dim_index] - subseq_size_ + 1)
            - synopsis_origin_[mbr_dim_index] + 1;
    cell_nums_[mbr_dim_index] = (mbr_dim_len - 1) / mbr_size_ + 1;
    /*
     * Turning off this check temporarily (?) since synopsis might cover only
     * a part of the original array. This is perfectly legal, but might get
     * the user into trouble, if she doesn't know this...
     */
//    if (mbr_dim.getCurrEnd() + 1 < cell_nums_[mbr_dim_index]) { // StartMin() == 0
//        std::ostringstream err_msg;
//        err_msg << "Synopsis must have at least "
//                << cell_nums_[mbr_dim_index] << " cells"
//                << " in the dimension " << mbr_dim.getBaseName();
//        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
//                SCIDB_LE_ILLEGAL_OPERATION) << err_msg.str();
//    }

    // Check if we have the required attributes
    if (!attributes_.count("low") || !attributes_.count("high")) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION) <<
                "DFT synopsis must have 'low' and 'high' attributes";
    }
}

IntervalValue SeqSynopsis::SqDist(const Coordinates &low,
        const Coordinates &high, const TransformedSequenceInfo &seq_info) {
    // Corrected coordinates
    Coordinates lowc{low}, highc{high};
    CheckAndCorrectBounds(lowc, highc);
    IntervalValue res; // NUL
    res.min_ = res.max_ = std::numeric_limits<double>::max();

    // Proper transformed sequence
    assert(seq_info.sequence_.count(subseq_size_));
    const auto &seq = seq_info.sequence_.find(subseq_size_)->second;
	assert(seq.size() % features_num_ == 0);

	// The sequence must fit
	if (highc.back() - lowc.back() + 1 < seq_info.original_length_) {
		// Query sequence doesn't fit -- return NULL
		return res;
	}
	/*
	 * Low/high should conform to the possible starting points.
	 */
	const size_t query_pieces = seq.size() / features_num_;
	highc.back() -= seq_info.original_length_ - query_pieces * subseq_size_;

	/*
	 * High actually defines the end of the search interval. Should adjust
	 * for the last subsequence, since we want low-high interval to determine
	 * possible starting positions.
	 */
	highc.back() -= subseq_size_ - 1;

	/*
	 * What we're going to do here is to go through all possible waveform
	 * combinations (assuming the last coordinates define the waveform
	 * interval, and the others determine possible different waveforms). For
	 * each waveform we measure the minimum sq. distance for each
	 * subsequence of the query sequence. Then, we sum them up and possibly
	 * apply the coefficient (e.g., in case of PAA). Note, since the
	 * individual distances are squared it is legal to sum them).
	 */
	struct SubSeqInfo {
	    RegionIterator iter_;
	    Coordinates low_, high_;
	    double min_sq_ = std::numeric_limits<double>::max();

	    SubSeqInfo(SeqSynopsis &syn, const Coordinates &low,
	            const Coordinates &high) :
	        iter_{syn, low, high}, low_{low}, high_{high} {}

	    // Reset based on the new _component_ coordinates (see the assert()).
        void Reset(const Coordinates &new_comp) {
	        assert(new_comp.size() + 1 == low_.size());
	        min_sq_ = std::numeric_limits<double>::max();
	        memcpy(low_.data(), new_comp.data(),
	            sizeof(Coordinate) * new_comp.size());
            memcpy(high_.data(), new_comp.data(),
                sizeof(Coordinate) * new_comp.size());
	        iter_.Reset(low_, high_);
	    }
	};
	std::vector<SubSeqInfo> subseqs;
	subseqs.reserve(query_pieces);

	// Iterator to go over components
	RowMajorIterator comp_it{Coordinates(lowc.begin(), lowc.end() - 1),
	    Coordinates(highc.begin(), highc.end() - 1)};
	// Init subsequence iterators
	Coordinates subseq_low{lowc}, subseq_high{lowc}; // the same component
	for (size_t i = 0; i < query_pieces; ++i) {
	    /*
	     *  Small optimization a certain subsequence might not start from
	     *  the beginning or go to the very end.
	     */
	    subseq_low.back() = lowc.back() + i * subseq_size_;
	    subseq_high.back() = highc.back() -
	            (query_pieces - i - 1) * subseq_size_;
	    subseqs.emplace_back(*this, subseq_low, subseq_high);
	}

	// Compute the result
    while (!comp_it.end()) {
        // Compute minsq for every subsequence
        bool comp_null = false;
        for (size_t i = 0; i < query_pieces; ++i) {
            auto &ss_info = subseqs[i];
            auto &it = ss_info.iter_;
            bool valid_cell = false;
            while (!it.end()) {
                const SeqCell &cell = it.GetCell();
                if (!cell.mbr_.low_.empty()) {
                    const double min_mbr_dist = cell.mbr_.MinSqDist(
                            seq.data() + i * features_num_);
                    if (min_mbr_dist < ss_info.min_sq_) {
                        ss_info.min_sq_ = min_mbr_dist;
                        valid_cell = true;
                    }
                }
                ++it;
            }
            if (!valid_cell) {
                // Only empty cells for the subsequence or NaN
                comp_null = true;
                break;
            }
        }
        // Here, we have computed every subsequence
        if (!comp_null) {
            res.state_ = IntervalValue::NON_NULL;
            double comp_sum = 0.0;
            for (size_t i = 0; i < query_pieces; ++i) {
                assert(subseqs[i].min_sq_ !=
                        std::numeric_limits<double>::max());
                comp_sum += subseqs[i].min_sq_;
            }
            if (comp_sum < res.min_) {
                res.min_ = comp_sum;
            }
        }
        // Move to the next component
        ++comp_it;
        if (!comp_it.end()) {
            // Reset point iterators
            for (auto &ss: subseqs) {
                ss.Reset(*comp_it);
            }
        }
    }
    if (res.state_ != IntervalValue::NUL && type_ == Type::PAA) {
        // PAA uses a little bit different lower-bound function
        res.min_ *= double(subseq_size_) / features_num_;
    }
	return res;
}

AggSynopsis::AggSynopsis(const ArrayDesc &data_desc,
                         const ArrayPtr &array) :
        Base{data_desc, array, array->getArrayDesc().getDimensions().size()} {
    /*
     * Additional check on the number of cells.
     */
    const ArrayDesc &synopsis_desc = array->getArrayDesc();
    const size_t dims = synopsis_desc.getDimensions().size();
    for (size_t i = 0; i < dims; ++i) {
        const DimensionDesc &syn_dim = synopsis_desc.getDimensions()[i];
        if (syn_dim.getCurrEnd() + 1 < cell_nums_[i]) { // StartMin() == 0
            std::ostringstream err_msg;
            err_msg << "Synopsis must have at least "
                    << cell_nums_[i] << " cells"
                    << " in the dimension " << syn_dim.getBaseName();
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                    SCIDB_LE_ILLEGAL_OPERATION) << err_msg.str();
        }
    }

    /*
     * Check if we have the required attributes.
     */
    const char * const attr_names[] = {"min", "max", "sum", "count"};
    for (const char *attr: attr_names) {
        const auto it = attributes_.find(attr);
        if (it == attributes_.end()) {
            std::ostringstream err_msg;
            err_msg << "Cannot find attribute " << attr << " in the sample "
                    << synopsis_array_->getArrayDesc().getName();
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                                   SCIDB_LE_ILLEGAL_OPERATION) << err_msg.str();
        }
    }
}

Sampler::Sampler(const ArrayDesc &data_desc,
        const std::vector<ArrayPtr> &synopsis_arrays,
        const SearchlightConfig &sl_config) :
            data_desc_{data_desc},
            sl_config_(sl_config),
            logger_(logger) {
    // Create catalog of describing synopses
    for (const auto &array: synopsis_arrays) {
        const std::string &array_name =
                ArrayDesc::makeUnversionedName(array->getName());
        const auto &params = TokenizeString(array_name, "_");

        // Attribute name is the second to last component of the name
        const std::string &attr_name = *(params.rbegin() + 1);
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
    cache_synopses_ = sl_config_.get("searchlight.sampler.cache_query", 1);
}

Sampler::~Sampler() {
    std::ostringstream stats_str;
    stats_str << "Stats for synopses follows:\n";
    for (size_t i = 0; i < synopses_.size(); i++) {
        stats_str << "\tStats for attribute " << i << ": \n";
        for (const auto &syn: synopses_[i].agg_synopses_) {
            syn->OutputStats(stats_str);
            stats_str << "\n";
        }
        for (const auto &syn: synopses_[i].seq_synopses_) {
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

void Sampler::LoadSampleForAttribute(const std::string &attr_name,
        AttributeID attr_search_id) {
    // Add a new vector of synopses if needed
    if (attr_search_id + 1 > synopses_.size()) {
        synopses_.resize(attr_search_id + 1);
    }

    // If we have loaded something for the attribute, ignore
    auto &attribute_synopses = synopses_[attr_search_id];
    if (attribute_synopses.loaded_) {
        return;
    }

    // See what we have for the attribute
    const auto &synops = array_synopses_[attr_name];
    for (const auto &syn: synops) {
    	// Determine if we're loading a DFT synopsis: they start as "dft_"
        const std::string &syn_name =
                ArrayDesc::makeUnversionedName(syn->getName());
        const std::string possible_type = TokenizeString(syn_name, "_")[0];
        if (possible_type == "dft" || possible_type == "paa") {
            const auto emplaced = seq_synopsis_cache_.emplace(syn_name,
                    data_desc_, syn);
            if (!emplaced.second) {
                // Found in the cache
                LOG4CXX_INFO(logger,
                    "Synopsis found in the cache: " << syn_name);
                emplaced.first->ReplaceArray(syn);
            }
            attribute_synopses.seq_synopses_.push_back(emplaced.first);
        } else {
            // Aggregate synopsis
            const auto emplaced = agg_synopsis_cache_.emplace(syn_name,
                    data_desc_, syn);
            if (!emplaced.second) {
                // Found in the cache
                LOG4CXX_INFO(logger,
                    "Synopsis found in the cache: " << syn_name);
                emplaced.first->ReplaceArray(syn);
            }
            attribute_synopses.agg_synopses_.push_back(emplaced.first);
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
    } else if (cache_type_string != "none") {
        LOG4CXX_WARN(logger, "Unknown cache type: " << cache_type_string);
    }

	// Are we preloading synopses?
    const bool preload_syns =
			sl_config_.get("searchlight.sampler.preload", 0);

    // Set cache type and preload aggregate synopses
    if (!attribute_synopses.agg_synopses_.empty()) {
		const size_t memory_limit_mb =
		        sl_config_.get("searchlight.sampler.memory_per_attr", 0UL);
    	PrepareSynopses(attribute_synopses.agg_synopses_, cache_type,
    	        preload_syns, memory_limit_mb, attr_name, attr_search_id);
    }

    // Set cache type and preload for DFT synopses
    if (!attribute_synopses.seq_synopses_.empty()) {
		const size_t memory_limit_mb =
				sl_config_.get("searchlight.sampler.memory_per_attr_seq", 0UL);
    	PrepareSynopses(attribute_synopses.seq_synopses_, cache_type,
    	        preload_syns, memory_limit_mb, attr_name, attr_search_id);
    }
}

bool AggSynopsis::RegionValidForSample(const Coordinates &low,
        const Coordinates &high) const {
    // copy to align
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

void AggSynopsis::ComputeAggregate(const Coordinates &low,
        const Coordinates &high, const SampleAggregatePtrVector &aggs) {
    // copy to check correctness and align
    Coordinates inlow(low), inhigh(high);
    CheckAndCorrectBounds(inlow, inhigh);

    // go through the chunks
    RegionIterator iter(*this, inlow, inhigh);
    uint64_t cells_accessed = 0;
    while (!iter.end()) {
        const AggCell &cell = iter.GetCell();
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

void AggSynopsis::ComputeAggregatesWithThr(
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
                const AggCell &cell = iter.GetCell();
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

IntervalValue AggSynopsis::GetElement(const Coordinates &point) {
    IntervalValue res;
    if (!CheckBounds(point)) {
        return res; // out of bounds -- NULL
    }

    RegionIterator iter(*this, point, point);
    const AggCell &cell = iter.GetCell();
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

    // Check if the region is valid (assume all synopses cover the same area)
    IntervalValueVector res(aggs.size());
    assert(!synopses_[s_attr].agg_synopses_.empty());
    if (!synopses_[s_attr].agg_synopses_[0]->CheckIfValid(low, high)) {
        // Not valid region, return nulls
        return res;
    }

    // Choose synopsis
    AggSynopsis *syn = nullptr;
    if (exact) {
        for (const auto &s: synopses_[s_attr].agg_synopses_) {
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
            // Clear res so that the caller can continue
            res.clear();
            return res;
        } else {
            syn->ComputeAggregate(low, high, aggs);
        }
    } else {
        // Non-exact computation: use the threshold algorithm
        assert(!synopses_[s_attr].agg_synopses_.empty());
        auto &syns = synopses_[s_attr].agg_synopses_;

        Region orig_region{low, high};
        syns[0]->CheckAndCorrectBounds(orig_region.low_, orig_region.high_);
        std::vector<Region> regions{orig_region};
        std::vector<Region> left_regions;

        for (size_t i = 0; i < syns.size(); i++) {
            // Check if we have a cached higher-resolution synopsis
            bool have_better_synopses = (i != syns.size() - 1) &&
                    syns[i + 1]->IsCached();

            // Check the MBR threshold
            if (mbr_thr_ > 0.0) {
                Region syn_mbr(orig_region);
                syns[i]->GetSynopsisMBR(syn_mbr.low_, syn_mbr.high_);
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
    for (size_t i = 0; i < aggs.size(); i++) {
        aggs[i].get()->Finalize(res[i]);
    }

    return res;
}

IntervalValue Sampler::GetElement(const Coordinates &point,
        AttributeID attr) const {
    // We always estimate element via the primary sample
    assert(!synopses_[attr].agg_synopses_.empty());
    const auto &syn = synopses_[attr].agg_synopses_[0];
    return syn->GetElement(point);
}

void Sampler::ComputeDFTs(const DoubleVector &seq, size_t ss_size,
		size_t dft_num, DoubleVector &res) {
    // Result space
    res.reserve(seq.size() / ss_size * dft_num);
	// In out arrays (reusable)
	double *in = (double *)fftw_malloc(sizeof(double) * ss_size);
	fftw_complex *out = (fftw_complex *)fftw_malloc(
			sizeof(fftw_complex) * (ss_size / 2 + 1));
	fftw_plan p = fftw_plan_dft_r2c_1d(ss_size, in, out, FFTW_ESTIMATE);
	// Break the sequence and compute DFTs
	const double sqrt_N = std::sqrt(ss_size);
	for (size_t i = 0; i + ss_size - 1 < seq.size(); i += ss_size) {
		// Copy to the planned array
		memcpy(in, seq.data() + i, sizeof(double) * ss_size);
		// Perform DFT
		fftw_execute(p);
		// Copy to the result
		for (size_t j = 0; j < dft_num / 2; ++j) {
			res.push_back(out[j][0] / sqrt_N);
			res.push_back(out[j][1] / sqrt_N);
		}
		if (dft_num % 2 != 0) {
			res.push_back(out[dft_num / 2][0] / sqrt_N);
		}
	}
	// Cleanup
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
}

void Sampler::ComputePAA(const DoubleVector &seq, size_t ss_size,
                         size_t feat_num, DoubleVector &res) {
    assert(ss_size % feat_num == 0);
    assert(seq.size() >= ss_size);
    const size_t win_size = ss_size / feat_num;
    res.reserve(seq.size() / ss_size * feat_num);
    for (size_t first = 0, last = first + ss_size - 1; last < seq.size();
            first += ss_size, last += ss_size) {
        for (size_t f = 0; f < feat_num; ++f) {
            double avg = 0.0;
            const size_t win_start = first + f * win_size;
            const size_t win_end = win_start + win_size - 1;
            for (size_t j = win_start; j <= win_end; ++j) {
                avg += seq[j];
            }
            avg /= win_size;
            res.push_back(avg);
        }
    }
}

void Sampler::RegisterQuerySequence(AttributeID attr, size_t seq_id,
		const DoubleVector &seq) {
    const auto ins = seq_cache_.emplace(seq_id,
            TransformedSequenceInfo{seq.size(), attr});
    assert(ins.second);
    auto &seq_info = ins.first->second;
	for (const auto &seq_syn: synopses_[attr].seq_synopses_) {
		const size_t ss_size = seq_syn->GetSubsequenceSize();
		const size_t feat_num = seq_syn->GetFeaturesNum();
		if (seq_syn->GetType() == SeqSynopsis::Type::PAA) {
		    ComputePAA(seq, ss_size, feat_num, seq_info.sequence_[ss_size]);
		} else {
		    ComputeDFTs(seq, ss_size, feat_num, seq_info.sequence_[ss_size]);
		}
	}
}

IntervalValue Sampler::SqDist(const Coordinates &low, const Coordinates &high,
        size_t seq_id) const {
    const auto it = seq_cache_.find(seq_id);
    assert(it != seq_cache_.end());
    const TransformedSequenceInfo &seq_info = it->second;
    // FIXME: For now we compute distance via the primary sample
    assert(!synopses_[seq_info.sattr_].seq_synopses_.empty());
    const auto &syn = synopses_[seq_info.sattr_].seq_synopses_[0];
    return syn->SqDist(low, high, seq_info);
}

void Sampler::ClearPersistentCache() {
    LOG4CXX_INFO(logger, "Clearing cache for Sampler...");
    agg_synopsis_cache_.Clear();
    seq_synopsis_cache_.Clear();
}

// Synopsis cache
Sampler::AggSynopsisCache Sampler::agg_synopsis_cache_;
Sampler::SeqSynopsisCache Sampler::seq_synopsis_cache_;

} /* namespace searchlight */

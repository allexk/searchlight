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
 * @file adapter.cpp
 * The implementation of the adapter.
 *
 * @author Alexander Kalinin
 */

#include "searchlight.h"
#include "adapter.h"
#include "array_desc.h"

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.adapter"));

Adapter::~Adapter() {
    std::ostringstream os;
    os << "Adapter(" << name_ << ") stats:\n";
    for (size_t i = 0; i < usage_stats_.size(); ++i) {
        const UsageStats &us = usage_stats_[i];
        const float secs =
                std::chrono::duration_cast<std::chrono::duration<double>>(
                us.total_req_time_).count();
        os << "Frame (" << i << ") time: " << secs << "s, accesses: D(" <<
                us.accesses_[3] << ") I(" << us.accesses_[2] << ") E(" <<
                us.accesses_[0] << ")\n";
    }
    LOG4CXX_INFO(logger, os.str());
}

void Adapter::UpdateStatsWithRegion(
        const Coordinates &low, const Coordinates &high) const {
    RegionIterator iter{array_desc_.GetOriginalArrayDesc(), low, high};
    while (!iter.end()) {
        stats_.chunks_pos_.insert(iter.CurrentPosition());
        ++iter;
    }
}

IntervalValueVector Adapter::ComputeAggregate(const Coordinates &low,
        const Coordinates &high, AttributeID attr,
        const StringVector &aggr_names) const {

    if (logger->isTraceEnabled()) {
        std::ostringstream deb_str;
        deb_str << "Aggregate request: Low: ";
        std::copy(low.begin(), low.end(),
                std::ostream_iterator<Coordinate>(deb_str, ", "));
        deb_str << " High: ";
        std::copy(high.begin(), high.end(),
                std::ostream_iterator<Coordinate>(deb_str, ", "));
        deb_str << " Attr: " << attr << " Aggs: ";
        std::copy(aggr_names.begin(), aggr_names.end(),
                std::ostream_iterator<std::string>(deb_str, ", "));
        deb_str << " Mode: " << mode_;

        logger->trace(deb_str.str(), LOG4CXX_LOCATION);
    }

    if (stats_enabled_) {
        UpdateStatsWithRegion(low, high);
    }

    // starting the timer
    const auto req_start_time = std::chrono::steady_clock::now();

    IntervalValueVector res; // NULLs by default
    if (mode_ == EXACT) {
        // First, try to use a synopsis
        const Sampler &sampler = array_desc_.GetSampler();
        res = sampler.ComputeAggregate(low, high, attr, aggr_names, true);
        if (res.empty()) {
            // Compute from real data
            const ArrayAccess &array = array_desc_.GetAccessor();
            TypedValueVector value_res = array.ComputeAggreagate(low, high,
                    attr, aggr_names);

            // Convert to the interval format
            res.resize(aggr_names.size()); // NULLs by default
            assert(value_res.size() == res.size());
            for (size_t i = 0; i < value_res.size(); i++) {
                const TypedValue &val_type = value_res[i];
                if (!val_type.first.isNull()) {
                   IntervalValue &r = res[i];
                   r.state_ = IntervalValue::NON_NULL;
                   r.min_ = r.max_ = r.val_ = scidb::ValueToDouble(
                           val_type.second.typeId(),
                           val_type.first);
                }
            }
        } else {
            LOG4CXX_TRACE(logger, "Computed from synopsis in the exact mode");
        }
    } else if (mode_ == APPROX || mode_ == INTERVAL) {
        const Sampler &sampler = array_desc_.GetSampler();
        res = sampler.ComputeAggregate(low, high, attr, aggr_names, false);
        assert(!res.empty());
        if (mode_ == APPROX) {
            for (IntervalValueVector::iterator it = res.begin();
                    it != res.end(); it++) {
                it->min_ = it->max_ = it->val_;
                if (it->state_ == IntervalValue::MAY_NULL) {
                    it->state_ = IntervalValue::NON_NULL;
                }
            }
        }
    } else if (mode_ == DUMB) {
        // We are being dumb: (-inf; +inf) and it may be NULL
        res.resize(aggr_names.size()); // NULLs by default
        for (auto &int_val: res) {
            int_val.min_ = std::numeric_limits<double>::lowest();
            int_val.max_ = std::numeric_limits<double>::max();
            int_val.state_ = IntervalValue::MAY_NULL;
            int_val.val_ = 0;
        }
    } else {
        // cannot happen
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Unknown adapter mode!";
    }

    // "stopping" the timer
    const auto req_end_time = std::chrono::steady_clock::now();
    usage_stats_.back().total_req_time_ +=
            std::chrono::duration_cast<decltype(UsageStats::total_req_time_)>(
            req_end_time - req_start_time);
    usage_stats_.back().accesses_[int(mode_)]++;

    if (logger->isTraceEnabled()) {
        std::ostringstream deb_str;
        deb_str << "Computed aggregates: ";
        std::copy(res.begin(), res.end(),
                std::ostream_iterator<IntervalValue>(deb_str, ", "));
        logger->trace(deb_str.str(), LOG4CXX_LOCATION);
    }

    return res;
}

IntervalValue Adapter::GetElement(const Coordinates &point,
        AttributeID attr) const {

    if (logger->isTraceEnabled()) {
        std::ostringstream deb_str;
        deb_str << "Element request: Point: ";
        std::copy(point.begin(), point.end(),
                std::ostream_iterator<Coordinate>(deb_str, ", "));
        deb_str << " Attr: " << attr;
        deb_str << " Mode: " << mode_;

        logger->trace(deb_str.str(), LOG4CXX_LOCATION);
    }

    if (stats_enabled_) {
        UpdateStatsWithRegion(point, point);
    }

    // starting the timer
    const auto req_start_time = std::chrono::steady_clock::now();

    IntervalValue res; // NULL
    if (mode_ == EXACT) {
        // First, resolve the attribute
        AttributeID real_id = array_desc_.GetArrayAttrributeID(attr);

        // Compute
        const ArrayAccess &array = array_desc_.GetAccessor();
        TypedValue tres;
        bool null = array.GetElement(point, real_id, tres);
        if (!null && !tres.first.isNull()) {
            res.state_ = IntervalValue::NON_NULL;
            res.min_ = res.max_ = res.val_ = scidb::ValueToDouble(
                tres.second.typeId(),
                tres.first);
        }
    } else if (mode_ == APPROX || mode_ == INTERVAL) {
        const Sampler &sampler = array_desc_.GetSampler();
        res = sampler.GetElement(point, attr);
        if (mode_ == APPROX) {
            res.min_ = res.max_ = res.val_;
            if (res.state_ == IntervalValue::MAY_NULL) {
                // TODO: make the choice with some probability returned in res?
                res.state_ = IntervalValue::NON_NULL;
            }
        }
    } else if (mode_ == DUMB) {
        // We are being dumb: (-inf; +inf) and it may be NULL
        res.min_ = std::numeric_limits<double>::lowest();
        res.max_ = std::numeric_limits<double>::max();
        res.state_ = IntervalValue::MAY_NULL;
        res.val_ = 0;
    } else {
        // cannot happen
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Unknown adapter mode!";
    }

    // "stopping" the timer
    const auto req_end_time = std::chrono::steady_clock::now();
    usage_stats_.back().total_req_time_ +=
            std::chrono::duration_cast<decltype(UsageStats::total_req_time_)>(
            req_end_time - req_start_time);
    usage_stats_.back().accesses_[int(mode_)]++;
    LOG4CXX_TRACE(logger, "Computed element: " << res);

    return res;
}

IntervalValue Adapter::SqDist(const Coordinates &low, const Coordinates &high,
		size_t int_dim,
		AttributeID attr, size_t seq_id) const {
	// Tracing
    if (logger->isTraceEnabled()) {
        std::ostringstream deb_str;
        deb_str << "Square distance: Low: ";
        std::copy(low.begin(), low.end(),
                std::ostream_iterator<Coordinate>(deb_str, ", "));
        deb_str << " High: ";
        std::copy(high.begin(), high.end(),
                std::ostream_iterator<Coordinate>(deb_str, ", "));
        deb_str << " Attr: " << attr;
        deb_str << " Seq ID: " << seq_id;
        deb_str << " Mode: " << mode_;
        logger->trace(deb_str.str(), LOG4CXX_LOCATION);
    }

    // Simulation stats
    if (stats_enabled_) {
        UpdateStatsWithRegion(low, high);
    }
    // starting the timer
    const auto req_start_time = std::chrono::steady_clock::now();

    IntervalValue res; // NULL
    if (mode_ == EXACT) {
        // First, resolve the attribute
        AttributeID real_attr_id = array_desc_.GetArrayAttrributeID(attr);
        // Compute
        const ArrayAccess &array = array_desc_.GetAccessor();
        double dist;
        const DoubleVector &query_seq = array_desc_.GetQuerySequence(seq_id);
        const bool success = array.SqDistance(low, high, int_dim, real_attr_id,
        		query_seq, dist);
        if (success) {
        	res.state_ = IntervalValue::NON_NULL;
        	res.min_ = res.max_ = res.val_ = dist;
        }
    } else if (mode_ == APPROX || mode_ == INTERVAL) {
        const Sampler &sampler = array_desc_.GetSampler();
        res = sampler.SqDist(attr, low[int_dim], high[int_dim], seq_id);
        if (mode_ == APPROX) {
            res.val_ = res.max_ = res.min_;
            if (res.state_ == IntervalValue::MAY_NULL) {
                // TODO: make the choice with some probability returned in res?
                res.state_ = IntervalValue::NON_NULL;
            }
        }
    } else if (mode_ == DUMB) {
        // We are being dumb: (0; +inf) and it may be NULL
        res.min_ = 0;
        res.max_ = std::numeric_limits<double>::max();
        res.state_ = IntervalValue::NON_NULL;
        res.val_ = 0;
    } else {
        // cannot happen
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Unknown adapter mode!";
    }

    // "stopping" the timer
    const auto req_end_time = std::chrono::steady_clock::now();
    usage_stats_.back().total_req_time_ +=
            std::chrono::duration_cast<decltype(UsageStats::total_req_time_)>(
            req_end_time - req_start_time);
    usage_stats_.back().accesses_[int(mode_)]++;

    LOG4CXX_TRACE(logger, "Computed square distance: " << res);

    return res;
}

void Adapter::SetCustomFail() const {
    if (sl_solver_id_ != -1) {
        sl_.GetSLSolver(sl_solver_id_).BeginCustomFail();
    }
}

std::ostream &operator<<(std::ostream &os, const IntervalValue &iv) {
    if (iv.state_ == IntervalValue::NUL) {
        os << "(NULL)";
    } else {
        os << "(Min: " << iv.min_ << ", Max: " << iv.max_ <<
                ", Val: " << iv.val_;
        if (iv.state_ == IntervalValue::MAY_NULL) {
            os << ", null?";
        }
        os << ")";
    }
    return os;
}
} /* namespace searchlight */

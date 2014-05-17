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

#include "adapter.h"

namespace searchlight {
IntervalValueVector Adapter::ComputeAggregate(const Coordinates &low,
        const Coordinates &high, AttributeID attr,
        const StringVector &aggr_names) const {

    IntervalValueVector res(aggr_names.size()); // NULLs by default
    if (mode_ == EXACT) {
        // First, resolve the attribute
        AttributeID real_id = array_desc_.GetArrayAttrributeID(attr);

        // Compute
        const ArrayAccess &array = array_desc_.GetAccessor();
        TypedValueVector value_res = array.ComputeAggreagate(low, high, real_id,
                aggr_names);

        // Convert to the interval format
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
    } else if (mode_ == APPROX || mode_ == INTERVAL) {
        const Sampler &sampler = array_desc_.GetSampler();
        res = sampler.ComputeAggregate(low, high, attr, aggr_names);
        if (mode_ == APPROX) {
            for (IntervalValueVector::iterator it = res.begin();
                    it != res.end(); it++) {
                it->min_ = it->max_ = it->val_;
                if (it->state_ == IntervalValue::MAY_NULL) {
                    it->state_ = IntervalValue::NON_NULL;
                }
            }
        }
    } else {
        // cannot happen
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Unknown adapter mode!";
    }

    return res;
}

IntervalValue Adapter::GetElement(const Coordinates &point,
        AttributeID attr) const {

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
                res.state_ = IntervalValue::NON_NULL;
            }
        }
    } else {
        // cannot happen
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Unknown adapter mode!";
    }

    return res;
}

} /* namespace searchlight */

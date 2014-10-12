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
 * @file searchlight_udfs.cpp
 * This file contains implementation for some useful standard "UDFs".
 *
 * @author Alexander Kalinin
 */

#include "ortools_inc.h"
#include "searchlight.h"

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.udfs"));

// A couple of utility functions
namespace {
    // Rounds a right bound up
    int64 RBoundToInt64(double rbound) {
        if (rbound == std::numeric_limits<double>::max()) {
            return std::numeric_limits<int64>::max();
        } else {
            return ceil(rbound);
        }
    }

    // Rounds a left bound down
    int64 LBoundToInt64(double lbound) {
        if (lbound == std::numeric_limits<double>::lowest()) {
            return std::numeric_limits<int64>::min();
        } else {
            return floor(lbound);
        }
    }
}

/**
 * This class allows to compute aggregate functions on arrays that are
 * expressed in terms of integer variables
 *
 * This is basically an extension of the IntExpr to properly handle
 * aggregated function over arrays.
 */
class AggrFuncExpr : public BaseIntExpr {
public:
    /**
     * The type of the aggregate.
     */
    enum AggType {
        SUM, //!< SUM the sum of elements within
        AVG, //!< AVG the average of elements within
        MIN, //!< MIN the minimum of elements within
        MAX, //!< MAX the maximum of elements within
    };

    /**
     * Constructs an aggregate expression for an array.
     *
     * @param s the solver we are using
     * @param agg the type of the aggregate
     * @param adapter the adapter for data access
     * @param low_lens_coords an array of coordinates + lengths for the region
     * @param params a vector of parameters; currently one -- the attribute
     */
    AggrFuncExpr(Solver *s, AggType agg, AdapterPtr adapter,
            const std::vector<IntVar *> &low_lens_coords,
            const std::vector<int64> &params) :
        BaseIntExpr(s),
        adapter_(adapter),
        attr_(AttributeID(params[0])),
        func_(agg),
        dims_(low_lens_coords.size() / 2), // assume even division (check later)
        low_lens_(low_lens_coords),
        low_lens_iters_(2 * dims_),
        min_(0),
        max_(0),
        min_support_low_(dims_),
        min_support_lens_(dims_),
        max_support_low_(dims_),
        max_support_lens_(dims_),
        min_max_init_(false) {

        if (low_lens_coords.size() != dims_ * 2 || params.size() != 1) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Aggregate UDFs coordinates/lengths are inconsistent "
                    "or parameters are incorrect!";
        }

        for (size_t i = 0; i < 2 * dims_; i++) {
            low_lens_iters_[i] = low_lens_[i]->MakeDomainIterator(true);
        }
    }

    /**
     * Destructor.
     */
    virtual ~AggrFuncExpr() {}

    /**
     * Returns the minimum of the expression among all windows.
     *
     * @return the minimum of the expression
     */
    virtual int64 Min() const {
        if (!ComputeMinMax()) {
            solver()->Fail();
        }
        return min_;
    }

    /**
     * Sets the minimum of the expression.
     *
     * Currently, we do not propagate it to the variables, but we fail the
     * search if the required minimum is impossible.
     *
     * @param m the minimum to set
     */
    virtual void SetMin(int64 m) {
        if (!ComputeMinMax() || m > max_) {
            solver()->Fail();
        }
    }

    /**
     * Returns the maximum of the expression among all windows.
     *
     * @return the maximum of the expression
     */
    virtual int64 Max() const {
        if (!ComputeMinMax()) {
            solver()->Fail();
        }
        return max_;
    }

    /**
     * Sets the maximum of the expression.
     *
     * Currently, we do not propagate it to the variables, but we fail the
     * search if the required maximum is impossible.
     *
     * @param m the maximum to set
     */
    virtual void SetMax(int64 m) {
        if (!ComputeMinMax() || m < min_) {
            solver()->Fail();
        }
    }

    /**
     * Attaches a demon to watch for range changes of the expression.
     *
     * @param d the demon to attach
     */
    virtual void WhenRange(Demon* d) {
        for (size_t i = 0; i < 2 * dims_; i++) {
            low_lens_[i]->WhenRange(d);
        }
    }

    /**
     * Returns a string representation for debug printing.
     *
     * @return a debug string
     */
    virtual std::string DebugString() const {
        std::string debug_str("AggregateArray(");
        for (int i = 0; i < dims_; i++) {
            if (i > 0) {
                debug_str += ", ";
            }
            debug_str += StringPrintf("%s + %s",
                    low_lens_[i]->DebugString().c_str(),
                    low_lens_[dims_ + i]->DebugString().c_str());
        }
        debug_str += ")";
        return debug_str;
    }

    /**
     * Accepts the given visitor for stats, buffering, printing, etc.
     *
     * @param visitor calling visitor
     */
    virtual void Accept(ModelVisitor* const visitor) const {
        std::string tag;
        switch (func_) {
            case SUM:
                tag = "UDF_sum";
                break;
            case AVG:
                tag = "UDF_avg";
                break;
            case MIN:
                tag = "UDF_min";
                break;
            case MAX:
                tag = "UDF_max";
                break;
        }

        visitor->BeginVisitIntegerExpression(tag, this);

        visitor->VisitIntegerVariableArrayArgument(ModelVisitor::kVarsArgument,
                low_lens_);

        std::vector<int64> params(1, int64(attr_));
        visitor->VisitIntegerArrayArgument(ModelVisitor::kValuesArgument,
                params);

        visitor->EndVisitIntegerExpression(tag, this);
    }

private:
    /*
     * Threshold of computing min/max values based on the individual regions
     * estimation instead of MBR-based.
     */
    static const int INDIVIDUAL_CHECK_THRESHOLD = 32;

    /*
     * Computes min/max of the aggregate if necessary.
     * Made const since it is called from Min()/Max(), which are made const
     * by or-tool's design.
     *
     * Returns true, if the min/max are set; false, if we need to backtrack.
     * We fail outside the function to avoid problems with setjmp/longjmp, since
     * (a) that's how Fail() is implemented (b) we have local objects there that
     * need dtors.
     */
    bool ComputeMinMax() const;

    /*
     *  Computes the aggregate for the given fixed coordinates.
     */
    IntervalValue ComputeFunc(const Coordinates &low, const Coordinates &high) const;

    /*
     * Estimates the aggregate for possible sub-regions of size
     * [min_size, max_size] of the specified region.
     */
    IntervalValue ComputeFuncSub(const Coordinates &low,
            const Coordinates &high, uint64_t min_size, uint64_t max_size) const;

    // Check if the support for min/max is valid.
    bool CheckSupport() const;

    // Returns the stringified aggregate name
    std::string GetAggregateName() const {
        std::string name;
        switch (func_) {
            case SUM:
                name = "sum";
                break;
            case AVG:
                name = "avg";
                break;
            case MIN:
                name = "min";
                break;
            case MAX:
                name = "max";
                break;
        }
        return name;
    }

    /*
     * Saves the old value at addr and assigns the new value
     * at the same place. Introduced to fix a small incompatibility
     * between int64 types in or-tools and scidb, which results in
     * ambiguity in SaveAndSetValue function.
     *
     * Strictly speaking, converting the pointer to (int64 *) is not
     * cool, but since they are both 8 bytes (guaranteed), it should be
     * fine.
     */
    void SaveCoordinate(Coordinate *addr, Coordinate new_val) const {
        solver()->SaveAndSetValue(reinterpret_cast<int64 *>(addr),
                int64(new_val));
    }

    const AdapterPtr adapter_; // the adapter for data access
    const AttributeID attr_; // attribute we are computing
    const AggType func_; // type of the aggregate
    const size_t dims_;  // dimensionality

    // coordinates for the window and the iterators (left point + lengths)
    const std::vector<IntVar *> low_lens_;
    std::vector<IntVarIterator *> low_lens_iters_;

    // Min/max aggregate values and caches for support
    mutable int64 min_;
    mutable int64 max_;
    mutable Coordinates min_support_low_, min_support_lens_;
    mutable Coordinates max_support_low_, max_support_lens_;
    mutable bool min_max_init_;
};

IntervalValue AggrFuncExpr::ComputeFunc(const Coordinates &low,
        const Coordinates &high) const {
    return adapter_->ComputeAggregate(low, high, attr_,
            StringVector(1, GetAggregateName()))[0];
}

IntervalValue AggrFuncExpr::ComputeFuncSub(const Coordinates &low,
        const Coordinates &high, uint64_t min_size, uint64_t max_size) const {
    StringVector agg_names(2);
    agg_names[0] = "min";
    agg_names[1] = "max";
    IntervalValueVector min_max = adapter_->ComputeAggregate(low, high, attr_,
            agg_names);

    IntervalValue res;
    if (min_max[0].state_ == IntervalValue::NON_NULL &&
            min_max[1].state_ == IntervalValue::NON_NULL) {
        res.state_ = IntervalValue::NON_NULL;
    } else if (min_max[0].state_ == IntervalValue::NUL &&
            min_max[1].state_ == IntervalValue::NUL) {
        return res; // will return NUL
    } else {
        res.state_ = IntervalValue::MAY_NULL;
    }

    switch(func_) {
        case AVG:
            res.min_ = min_max[0].min_;
            res.max_ = min_max[1].max_;
            break;
        case SUM:
            res.min_ = min_max[0].min_ * min_size;
            res.max_ = min_max[1].max_ * max_size;
            break;
        case MIN:
            res.min_ = min_max[0].min_;
            res.max_ = min_max[1].max_;
            break;
        case MAX:
            res.min_ = min_max[0].min_;
            res.max_ = min_max[1].max_;
            break;
    }

    LOG4CXX_TRACE(logger, "Computed MBR value: " << res);

    return res;
}

bool AggrFuncExpr::CheckSupport() const {
    // first time: no support
    if (!min_max_init_) {
        return false;
    }

    // the support is not valid if it is not a valid window anymore
    for (size_t i = 0; i < dims_; i++) {
        if (!low_lens_[i]->Contains(min_support_low_[i]) ||
                !low_lens_[dims_ + i]->Contains(min_support_lens_[i]) ||
                !low_lens_[i]->Contains(max_support_low_[i]) ||
                !low_lens_[dims_ + i]->Contains(max_support_lens_[i])) {
            return false;
        }
    }
    return true;
}

bool AggrFuncExpr::ComputeMinMax() const {
    IntervalValue new_min_max;
    Coordinates new_min_support_low, new_min_support_lens;
    Coordinates new_max_support_low, new_max_support_lens;

    // First case: variables are bound
    bool vars_bound = true;
    for (size_t i = 0; i < 2 * dims_; i++) {
        if (!low_lens_[i]->Bound()) {
            vars_bound = false;
            break;
        }
    }

    if (vars_bound) {
        if (CheckSupport()) {
            return true;
        }
        Coordinates low(dims_), high(dims_), lens(dims_);
        for (size_t i = 0; i < dims_; i++) {
            low[i] = low_lens_[i]->Value();
            lens[i] = low_lens_[dims_ + i]->Value();
            high[i] = low[i] + lens[i] - 1;
        }

        new_min_max = ComputeFunc(low, high);
        new_min_support_low = new_max_support_low = low;
        new_min_support_lens = new_max_support_lens = lens;
    } else {
        // Second case: below threshold and individual checks
        uint64_t reg_num = 1;
        for (size_t i = 0; i < 2 * dims_; i++) {
            reg_num *= low_lens_[i]->Size();
        }

        if (reg_num <= INDIVIDUAL_CHECK_THRESHOLD) {
            if (CheckSupport()) {
                return true;
            }
            /*
             * Here we need to go through every possible combination
             * of coords/lens, which might be tricky for an arbitrary
             * number of dimensions...
             */
            Coordinates low(dims_), high(dims_), lens(dims_);
            for (size_t i = 0; i < dims_; i++) {
                low_lens_iters_[i]->Init();
                low_lens_iters_[dims_ + i]->Init();
                low[i] = low_lens_iters_[i]->Value();
                lens[i] = low_lens_iters_[dims_ + i]->Value();
                high[i] = low[i] + lens[i] - 1;
            }

            while (true) {
                const IntervalValue val = ComputeFunc(low, high);
                if (new_min_max.state_ == IntervalValue::NUL) {
                    new_min_max = val;
                    new_min_support_low = new_max_support_low = low;
                    new_min_support_lens = new_max_support_lens = high;
                } else if (val.state_ != IntervalValue::NUL) {
                    if (val.min_ < new_min_max.min_) {
                        new_min_max.min_ = val.min_;
                        new_min_support_low = low;
                        new_min_support_lens = high;
                    }
                    if (val.max_ > new_min_max.max_) {
                        new_min_max.max_ = val.max_;
                        new_max_support_low = low;
                        new_max_support_lens = high;
                    }
                    if (val.state_ == IntervalValue::MAY_NULL) {
                        new_min_max.state_ = IntervalValue::MAY_NULL;
                    }
                }

                /*
                 *  Move to a new region. We try to move the current
                 *  iterator while possible, and then reset and move to
                 *  the next one and so on...
                 */
                size_t i = 0;
                while (i < 2 * dims_) {
                    IntVarIterator *it = low_lens_iters_[i];
                    const size_t k = i;

                    it->Next();
                    if (!it->Ok()) {
                        it->Init();
                        i++;
                    }

                    if (k < dims_) {
                        low[k] = it->Value();
                        high[k] = low[k] + lens[k] - 1;
                    } else {
                        lens[k - dims_] = it->Value();
                        high[k - dims_] = lens[k - dims_] +
                                low[k - dims_] - 1;
                    }

                    if (k == i) {
                        break;
                    }
                }

                // more regions? not if we have exhausted all iterators
                if (i == 2 * dims_) {
                    break;
                }
            }
        } else {
            // Third case: MBR + min/max based
            Coordinates low(dims_), high(dims_), lens(dims_);
            uint64_t min_size = 1, max_size = 1;
            bool mbr_changed = false;
            for (size_t i = 0; i < dims_; i++) {
                low[i] = low_lens_[i]->Min();
                high[i] = low_lens_[i]->Max() +
                        low_lens_[dims_ + i]->Max() - 1;
                lens[i] = high[i] - low[i] + 1;
                min_size *= low_lens_[dims_ + i]->Min();
                max_size *= low_lens_[dims_ + i]->Max();

                // min/max supports are equal for MBRs
                if (low[i] != min_support_low_[i] ||
                        lens[i] != min_support_lens_[i]) {
                    mbr_changed = true;
                }
            }

            if (min_max_init_ && !mbr_changed) {
                /*
                 * Strictly speaking, we might have different
                 * min_size or max_size here, but this is probably
                 * going to be really rare.
                 */
                return true;
            }

            new_min_max = ComputeFuncSub(low, high, min_size, max_size);
            new_min_support_low = new_max_support_low = low;
            new_min_support_lens = new_max_support_lens = lens;
        }
    }

    /*
     * At this point we have a proper new_min_max and supports. We need
     * to analyze the value and properly save it.
     */
    Solver * const s = solver();
    if (new_min_max.state_ == IntervalValue::NUL) {
        return false;
    }

    /*
     * Need to convert from double to int64, since or-tools do not support
     * floating values. For now, just round them to the nearest integer for the
     * exact value case or "nearest" interval for the interval case.
     *
     * TODO: revisit later if floats/doubles become available in or-tools.
     */
    int64 new_min = LBoundToInt64(new_min_max.min_);
    int64 new_max = RBoundToInt64(new_min_max.max_);
    if (new_min_max.min_ == new_min_max.max_) {
        new_min = new_max = round(new_min_max.min_);
    }

    // save the values and supports
    s->SaveAndSetValue(&min_, new_min);
    s->SaveAndSetValue(&max_, new_max);
    s->SaveAndSetValue(&min_max_init_, true);
    for (size_t i = 0; i < dims_; i++) {
        SaveCoordinate(&min_support_low_[i], new_min_support_low[i]);
        SaveCoordinate(&min_support_lens_[i], new_min_support_lens[i]);
        SaveCoordinate(&max_support_low_[i], new_max_support_low[i]);
        SaveCoordinate(&max_support_lens_[i], new_max_support_lens[i]);
    }

    return true;
}

extern "C"
IntExpr *CreateUDF_avg(Solver *solver, AdapterPtr adapter,
        const std::vector<IntVar *> &coord_lens,
        const std::vector<int64> &params) {
    return new AggrFuncExpr(solver, AggrFuncExpr::AVG, adapter, coord_lens,
            params);
}

extern "C"
IntExpr *CreateUDF_sum(Solver *solver, AdapterPtr adapter,
        const std::vector<IntVar *> &coord_lens,
        const std::vector<int64> &params) {
    return new AggrFuncExpr(solver, AggrFuncExpr::SUM, adapter, coord_lens,
            params);
}

extern "C"
IntExpr *CreateUDF_min(Solver *solver, AdapterPtr adapter,
        const std::vector<IntVar *> &coord_lens,
        const std::vector<int64> &params) {
    return new AggrFuncExpr(solver, AggrFuncExpr::MIN, adapter, coord_lens,
            params);
}

extern "C"
IntExpr *CreateUDF_max(Solver *solver, AdapterPtr adapter,
        const std::vector<IntVar *> &coord_lens,
        const std::vector<int64> &params) {
    return new AggrFuncExpr(solver, AggrFuncExpr::MAX, adapter, coord_lens,
            params);
}

} /* namespace searchlight */

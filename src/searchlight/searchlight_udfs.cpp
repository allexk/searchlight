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

#include "searchlight_udfs.h"
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

SearchlightUDF::SearchlightUDF(Solver *s, State *state,
        const AdapterPtr &adapter, AttributeID attr,
        const std::vector<IntExpr *> &vars) :
                BaseIntExpr{s},
                adapter_{adapter},
                attr_{attr},
                udf_state_{state},
                orig_vars_{vars} {
    // Init parameter vars
    param_vars_.reserve(vars.size());
    for (auto v: vars) {
        param_vars_.emplace_back(v);
    }
}

SearchlightUDF::ParameterVar::ParameterVar(const IntExpr *expr) :
        real_var_(dynamic_cast<const IntVar *>(expr) != nullptr),
        expr_{expr} {
    // Init iterator
    if (real_var_) {
        iter_ = static_cast<const IntVar *>(expr_)->MakeDomainIterator(true);
    } else {
        iter_ = expr->solver()->RevAlloc(new RangeIterator(this));
    }
}

bool SearchlightUDF::ParameterVar::Contains(int64 x) const {
    if (x < min_ || x > max_) {
        return false;
    }
    if (real_var_) {
        return static_cast<const IntVar *>(expr_)->Contains(x);
    } else {
        return x >= min_ && x <= max_;
    }
}

uint64 SearchlightUDF::ParameterVar::Size() const {
    if (real_var_) {
        return static_cast<const IntVar *>(expr_)->Size();
    } else {
        return max_ - min_ + 1;
    }
}

void SearchlightUDF::ParameterVar::Update() {
    min_ = expr_->Min();
    max_ = expr_->Max();
}

void SearchlightUDF::VisitOriginalExprs(ModelVisitor *visitor) const {
    visitor->VisitIntegerArgument("expr_size", int64(orig_vars_.size()));
    for (size_t i = 0; i < orig_vars_.size(); ++i) {
        const std::string etag = "expr_" + std::to_string(i);
        visitor->VisitIntegerExpressionArgument(etag, orig_vars_[i]);
    }
}

/**
 * Square distance from the array values to the query sequence.
 *
 * This is basically an extension of the IntExpr to properly handle
 * distance function over arrays.
 */
class SqDistFuncExpr : public SearchlightUDF {
public:
    /**
     * Constructs a distance function expression for an array.
     *
     * @param s the solver we are using
     * @param adapter the adapter for data access
     * @param coords array coordinates
     * @param params array: attribute, seq. id, seq. length, seq. dimension
     */
	SqDistFuncExpr(Solver *s, AdapterPtr adapter,
            const std::vector<IntExpr *> &coords,
            const std::vector<int64> &params) :
        SearchlightUDF(s, &state_, adapter, AttributeID(params[0]), coords),
        dims_(coords.size()),
		seq_id_(size_t(params[1])),
		seq_len_(size_t(params[2])),
		state_(dims_) {}

    /**
     * Returns a string representation for debug printing.
     *
     * @return a debug string
     */
    virtual std::string DebugString() const override {
        const auto &vars = GetVars();
        std::string debug_str("SqDistArray(");
        for (int i = 0; i < dims_; i++) {
            if (i > 0) {
                debug_str += ", ";
            }
            debug_str += StringPrintf("%s",
                    vars[i].DebugString().c_str());
        }
        debug_str += ")";
        return debug_str;
    }

    /**
     * Accepts the given visitor for stats, buffering, printing, etc.
     *
     * @param visitor calling visitor
     */
    virtual void Accept(ModelVisitor* const visitor) const override {
    	// Function start
        std::string tag("UDF_sqdist");
        visitor->BeginVisitIntegerExpression(tag, this);
        // Variables
        VisitOriginalExprs(visitor);
        // Params
        std::vector<int64> prms{int64(attr_), int64(seq_id_), int64(seq_len_)};
        visitor->VisitIntegerArrayArgument(ModelVisitor::kValuesArgument, prms);
        // Function end
        visitor->EndVisitIntegerExpression(tag, this);
    }

    /**
     * Saves state.
     */
    virtual const SearchlightUDF::State *SaveState() const override {
        if (state_.min_max_init_) {
            const SearchlightUDF::State *saved_state = new SqDistState(state_);
            return saved_state;
        } else {
            return nullptr;
        }
    }

    /**
     * Load state from the buf.
     *
     * @param state saved state
     */
    virtual void LoadState(const SearchlightUDF::State *state) override {
        if (state) {
            // Static cast is fine, it's controlled internal usage
            const SqDistState *saved_state =
                    static_cast<const SqDistState *>(state);
            Solver * const s = solver();
            s->SaveAndSetValue(&state_.min_, saved_state->min_);
            s->SaveAndSetValue(&state_.min_max_init_, true);
            for (size_t i = 0; i < dims_; i++) {
                SaveCoordinate(&state_.min_support_low_[i],
                        saved_state->min_support_low_[i]);
                SaveCoordinate(&state_.min_support_high_[i],
                    saved_state->min_support_high_[i]);
            }
        }
    }

private:
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
    virtual bool ComputeMinMax() const override;

    /*
     *  Computes the aggregate for the given fixed coordinates.
     */
    IntervalValue ComputeFunc(Coordinates &low, Coordinates &high) const;

    // Check if the support for min/max is valid.
    bool CheckSupport(const std::vector<ParameterVar> &vars) const;

    // Returns the stringified aggregate name
    std::string GetAggregateName() const {
    	return "UDF_sqdist";
    }

    const size_t dims_;  // dimensionality
    const size_t seq_id_; // Sequence ID
    const size_t seq_len_; // Sequence length

    // Min/max aggregate values and caches for support
    struct SqDistState : public SearchlightUDF::State {
        // Support for minimum: min, max values for vars
        Coordinates min_support_low_, min_support_high_;

        SqDistState(size_t dims) :
            min_support_low_(dims),
            min_support_high_(dims) {
            // Set max to a different value (its always max for this function)
            max_ = std::numeric_limits<int64>::max();
        }
    };
    mutable SqDistState state_;
};

IntervalValue SqDistFuncExpr::ComputeFunc(Coordinates &low,
        Coordinates &high) const {
    // SqDist requires high to determine the high point of the interval
    high.back() += seq_len_ - 1;
    const auto res = adapter_->SqDist(low, high, attr_, seq_id_);
    high.back() -= seq_len_ - 1;
    return res;
}

bool SqDistFuncExpr::CheckSupport(const std::vector<ParameterVar> &vars) const {
    // first time: no support
    if (!state_.min_max_init_) {
        return false;
    }

    /*
     * The support is valid if the variable intervals are the same.
     *
     * Strictly speaking this is not true, since the domains might have changed.
     * For example, the heuristic might have removed a value from the middle.
     * Min value is still value in this case, though, but possibly not as tight
     * as it could be.
     */
    for (size_t i = 0; i < dims_; i++) {
    	if (Coordinate(vars[i].Min()) != state_.min_support_low_[i] ||
    			Coordinate(vars[i].Max()) != state_.min_support_high_[i]) {
    		return false;
    	}
    }
    return true;
}

bool SqDistFuncExpr::ComputeMinMax() const {
    IntervalValue new_min;
    Coordinates new_min_support_low, new_min_support_high;
    const auto &vars = GetVars();

    // First case: variables are bound
    bool vars_bound = true;
    for (size_t i = 0; i < dims_; i++) {
        if (!vars[i].Bound()) {
            vars_bound = false;
            break;
        }
    }
    if (vars_bound) {
        if (CheckSupport(vars)) {
            return true;
        }
        Coordinates low(dims_), high(dims_);
        for (size_t i = 0; i < dims_; i++) {
            high[i] = low[i] = vars[i].Value();
        }
        // Compute the function
        new_min = ComputeFunc(low, high);
        new_min_support_low = low;
        new_min_support_high = high;
    } else {
        // Second case: checks for all combinations of non-sequence dimensions
		if (CheckSupport(vars)) {
			return true;
		}
		/*
		 * Here we need to go through every possible combination
		 * of non-sequence vars, which might be tricky for an arbitrary
		 * number of dimensions...
		 */
		Coordinates low(dims_), high(dims_);
		for (size_t i = 0; i < dims_ - 1; i++) {
            auto iter = vars[i].Iterator();
            iter->Init();
            high[i] = low[i] = iter->Value();
		}
		// For the sequence dimension keep it fixed at min-max
		low.back() = vars.back().Min();
		high.back() = vars.back().Max();

		while (true) {
			const IntervalValue val = ComputeFunc(low, high);
			if (new_min.state_ == IntervalValue::NUL) {
				new_min = val;
				new_min_support_low = low;
				new_min_support_high = high;
			} else {
				if (val.min_ < new_min.min_) {
					new_min.min_ = val.min_;
					new_min_support_low = low;
					new_min_support_high = high;
				}
			}

			/*
			 *  Move to a new region. We try to move the current
			 *  iterator while possible, and then reset and move to
			 *  the next one and so on... The sequence dimension is fixed.
			 */
			size_t i = 0;
			while (i < dims_ - 1) {
				// Try to move to the next value
				IntVarIterator *it = vars[i].Iterator();
				const size_t k = i;
				it->Next();
				if (!it->Ok()) {
					it->Init();
					i++;
				}
				high[k] = low[k] = it->Value();
				if (k == i) {
					// New combination: break and compute it
					break;
				}
			}
			// more regions? not if we have exhausted all iterators
			if (i == dims_ - 1) {
				break;
			}
		}
    }

    /*
     * At this point we have a proper new_min and support. We need
     * to analyze the value and properly save it.
     */
    Solver * const s = solver();
    int64 new_min_rounded;
    if (new_min.state_ == IntervalValue::NUL) {
        // If NUL, set the "never possible" value
        //new_min_rounded = kint64max;
        s->SaveAndSetValue(&state_.min_, kint64max);
        s->SaveAndSetValue(&state_.min_max_init_, false);
        return false;
    } else {
        /*
         * Need to convert from double to int64, since or-tools do not support
         * floating values. For now, just round them to the nearest integer for the
         * exact value case or "nearest" interval for the interval case.
         *
         * TODO: revisit later if floats/doubles become available in or-tools.
         */
        new_min_rounded = RBoundToInt64(new_min.min_);
    }

    // save the values and supports
    const bool can_cache = adapter_->CanCacheResults();
    s->SaveAndSetValue(&state_.min_, new_min_rounded);
    s->SaveAndSetValue(&state_.min_max_init_, can_cache);
    if (can_cache) {
        for (size_t i = 0; i < dims_; i++) {
            SaveCoordinate(&state_.min_support_low_[i], new_min_support_low[i]);
            SaveCoordinate(&state_.min_support_high_[i], new_min_support_high[i]);
        }
    }

    return true;
}

/**
 * This class allows to compute aggregate functions on arrays that are
 * expressed in terms of integer variables
 *
 * This is basically an extension of the IntExpr to properly handle
 * aggregated function over arrays.
 */
class AggrFuncExpr : public SearchlightUDF {
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
            const std::vector<IntExpr *> &low_lens_coords,
            const std::vector<int64> &params) :
        SearchlightUDF(s, &state_, adapter, AttributeID(params[0]), low_lens_coords),
        func_(agg),
        dims_(low_lens_coords.size() / 2), // assume even division (check later)
        state_(dims_) {

        if (low_lens_coords.size() != dims_ * 2 || params.size() != 1) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Aggregate UDFs coordinates/lengths are inconsistent "
                    "or parameters are incorrect!";
        }
    }

    /**
     * Destructor.
     */
    virtual ~AggrFuncExpr() {}

    /**
     * Returns a string representation for debug printing.
     *
     * @return a debug string
     */
    virtual std::string DebugString() const override {
        const auto &vars = GetVars();
        std::string debug_str("AggregateArray(");
        for (int i = 0; i < dims_; i++) {
            if (i > 0) {
                debug_str += ", ";
            }
            debug_str += StringPrintf("%s + %s",
                    vars[i].DebugString().c_str(),
                    vars[dims_ + i].DebugString().c_str());
        }
        debug_str += ")";
        return debug_str;
    }

    /**
     * Accepts the given visitor for stats, buffering, printing, etc.
     *
     * @param visitor calling visitor
     */
    virtual void Accept(ModelVisitor* const visitor) const override {
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
        // Expressions
        VisitOriginalExprs(visitor);
        // Parameters
        std::vector<int64> params(1, int64(attr_));
        visitor->VisitIntegerArrayArgument(ModelVisitor::kValuesArgument,
                params);
        visitor->EndVisitIntegerExpression(tag, this);
    }

    /**
     * Saves state.
     */
    virtual const SearchlightUDF::State *SaveState() const override {
        if (state_.min_max_init_) {
            const SearchlightUDF::State *saved_state = new AggState(state_);
            return saved_state;
        } else {
            return nullptr;
        }
    }

    /**
     * Load state from the state.
     *
     * @param state previously saved state
     */
    virtual void LoadState(const SearchlightUDF::State *state) override {
        if (state) {
            // Static cast is fine, it's controlled internal usage
            const AggState *saved_state = static_cast<const AggState *>(state);
            Solver * const s = solver();
            s->SaveAndSetValue(&state_.min_, saved_state->min_);
            s->SaveAndSetValue(&state_.max_, saved_state->max_);
            s->SaveAndSetValue(&state_.min_max_init_, true);
            for (size_t i = 0; i < dims_; i++) {
                SaveCoordinate(&state_.min_support_low_[i],
                        saved_state->min_support_low_[i]);
                SaveCoordinate(&state_.min_support_lens_[i],
                    saved_state->min_support_lens_[i]);
                SaveCoordinate(&state_.max_support_low_[i],
                    saved_state->max_support_low_[i]);
                SaveCoordinate(&state_.max_support_lens_[i],
                    saved_state->max_support_lens_[i]);
            }
        }
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
    virtual bool ComputeMinMax() const override;

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
    bool CheckSupport(const std::vector<ParameterVar> &vars) const;

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

    const AggType func_; // type of the aggregate
    const size_t dims_;  // dimensionality

    // Min/max aggregate values and caches for support
    struct AggState : public SearchlightUDF::State {
        Coordinates min_support_low_, min_support_lens_;
        Coordinates max_support_low_, max_support_lens_;

        AggState(size_t dims) :
            min_support_low_(dims),
            min_support_lens_(dims),
            max_support_low_(dims),
            max_support_lens_(dims) {}
    };
    mutable AggState state_;
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

bool AggrFuncExpr::CheckSupport(const std::vector<ParameterVar> &vars) const {
    // first time: no support
    if (!state_.min_max_init_) {
        return false;
    }

    // the support is not valid if it is not a valid window anymore
    for (size_t i = 0; i < dims_; i++) {
        if (!vars[i].Contains(state_.min_support_low_[i]) ||
                !vars[dims_ + i].Contains(state_.min_support_lens_[i]) ||
                !vars[i].Contains(state_.max_support_low_[i]) ||
                !vars[dims_ + i].Contains(state_.max_support_lens_[i])) {
            return false;
        }
    }
    return true;
}

bool AggrFuncExpr::ComputeMinMax() const {
    IntervalValue new_min_max;
    Coordinates new_min_support_low, new_min_support_lens;
    Coordinates new_max_support_low, new_max_support_lens;
    const auto &vars = GetVars();

    // First case: variables are bound
    bool vars_bound = true;
    for (size_t i = 0; i < 2 * dims_; i++) {
        if (!vars[i].Bound()) {
            vars_bound = false;
            break;
        }
    }

    if (vars_bound) {
        if (CheckSupport(vars)) {
            return true;
        }
        Coordinates low(dims_), high(dims_), lens(dims_);
        for (size_t i = 0; i < dims_; i++) {
            low[i] = vars[i].Value();
            lens[i] = vars[dims_ + i].Value();
            high[i] = low[i] + lens[i] - 1;
        }

        new_min_max = ComputeFunc(low, high);
        new_min_support_low = new_max_support_low = low;
        new_min_support_lens = new_max_support_lens = lens;
    } else {
        // Second case: below threshold and individual checks
        uint64_t reg_num = 1;
        for (size_t i = 0; i < 2 * dims_; i++) {
            reg_num *= vars[i].Size();
        }

        if (reg_num <= INDIVIDUAL_CHECK_THRESHOLD) {
            if (CheckSupport(vars)) {
                return true;
            }
            /*
             * Here we need to go through every possible combination
             * of coords/lens, which might be tricky for an arbitrary
             * number of dimensions...
             */
            Coordinates low(dims_), high(dims_), lens(dims_);
            for (size_t i = 0; i < dims_; i++) {
                vars[i].Iterator()->Init();
                vars[dims_ + i].Iterator()->Init();
                low[i] = vars[i].Iterator()->Value();
                lens[i] = vars[dims_ + i].Iterator()->Value();
                high[i] = low[i] + lens[i] - 1;
            }

            while (true) {
                const IntervalValue val = ComputeFunc(low, high);
                if (new_min_max.state_ == IntervalValue::NUL) {
                    new_min_max = val;
                    new_min_support_low = new_max_support_low = low;
                    new_min_support_lens = new_max_support_lens = lens;
                } else if (val.state_ != IntervalValue::NUL) {
                    if (val.min_ < new_min_max.min_) {
                        new_min_max.min_ = val.min_;
                        new_min_support_low = low;
                        new_min_support_lens = lens;
                    }
                    if (val.max_ > new_min_max.max_) {
                        new_min_max.max_ = val.max_;
                        new_max_support_low = low;
                        new_max_support_lens = lens;
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
                    IntVarIterator *it = vars[i].Iterator();
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
                    	// Successful move: i didn't increment
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
                low[i] = vars[i].Min();
                high[i] = vars[i].Max() +
                        vars[dims_ + i].Max() - 1;
                lens[i] = high[i] - low[i] + 1;
                min_size *= vars[dims_ + i].Min();
                max_size *= vars[dims_ + i].Max();

                // min/max supports are equal for MBRs
                if (low[i] != state_.min_support_low_[i] ||
                        lens[i] != state_.min_support_lens_[i]) {
                    mbr_changed = true;
                }
            }

            if (state_.min_max_init_ && !mbr_changed) {
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
    int64 new_min, new_max;
    if (new_min_max.state_ == IntervalValue::NUL) {
        // If NULL, set the "always fail" values
        //new_min = kint64max;
        //new_max = kint64min;
        s->SaveAndSetValue(&state_.min_, kint64max);
        s->SaveAndSetValue(&state_.max_, kint64min);
        s->SaveAndSetValue(&state_.min_max_init_, false);
        return false;
    } else {
        /*
         * Need to convert from double to int64, since or-tools do not support
         * floating values. For now, just round them to the nearest integer for the
         * exact value case or "nearest" interval for the interval case.
         *
         * TODO: revisit later if floats/doubles become available in or-tools.
         */
        if (new_min_max.min_ == new_min_max.max_) {
            new_min = new_max = round(new_min_max.min_);
        } else {
            new_min = LBoundToInt64(new_min_max.min_);
            new_max = RBoundToInt64(new_min_max.max_);
        }
    }

    // save the values and supports
    const bool can_cache = adapter_->CanCacheResults();
    s->SaveAndSetValue(&state_.min_, new_min);
    s->SaveAndSetValue(&state_.max_, new_max);
    s->SaveAndSetValue(&state_.min_max_init_, can_cache);
    if (can_cache) {
        for (size_t i = 0; i < dims_; i++) {
            SaveCoordinate(&state_.min_support_low_[i], new_min_support_low[i]);
            SaveCoordinate(&state_.min_support_lens_[i], new_min_support_lens[i]);
            SaveCoordinate(&state_.max_support_low_[i], new_max_support_low[i]);
            SaveCoordinate(&state_.max_support_lens_[i], new_max_support_lens[i]);
        }
    }

    return true;
}

extern "C"
IntExpr *CreateUDF_avg(Solver *solver, AdapterPtr adapter,
        const std::vector<IntExpr *> &coord_lens,
        const std::vector<int64> &params) {
    return new AggrFuncExpr(solver, AggrFuncExpr::AVG, adapter, coord_lens,
            params);
}

extern "C"
IntExpr *CreateUDF_sum(Solver *solver, AdapterPtr adapter,
        const std::vector<IntExpr *> &coord_lens,
        const std::vector<int64> &params) {
    return new AggrFuncExpr(solver, AggrFuncExpr::SUM, adapter, coord_lens,
            params);
}

extern "C"
IntExpr *CreateUDF_min(Solver *solver, AdapterPtr adapter,
        const std::vector<IntExpr *> &coord_lens,
        const std::vector<int64> &params) {
    return new AggrFuncExpr(solver, AggrFuncExpr::MIN, adapter, coord_lens,
            params);
}

extern "C"
IntExpr *CreateUDF_max(Solver *solver, AdapterPtr adapter,
        const std::vector<IntExpr *> &coord_lens,
        const std::vector<int64> &params) {
    return new AggrFuncExpr(solver, AggrFuncExpr::MAX, adapter, coord_lens,
            params);
}

extern "C"
IntExpr *CreateUDF_sqdist(Solver *solver, AdapterPtr adapter,
        const std::vector<IntExpr *> &coords,
        const std::vector<int64> &params) {
    return new SqDistFuncExpr(solver, adapter, coords, params);
}
} /* namespace searchlight */

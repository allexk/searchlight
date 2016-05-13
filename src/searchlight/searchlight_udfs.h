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
 * @file searchlight_udfs.h
 * Common definintions for standard UDFs.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_SEARCHLIGHT_UDFS_H_
#define SEARCHLIGHT_SEARCHLIGHT_UDFS_H_

#include "ortools_inc.h"
#include "scidb_inc.h"
#include "adapter.h"

namespace searchlight {

/**
 * Base class for all Searchlight UDFs.
 */
class SearchlightUDF : public BaseIntExpr {
public:
    /**
     * Abstract class for UDF states.
     */
    struct State {
        int64 min_ = 0;
        int64 max_ = 0;
        bool min_max_init_ = false;

        virtual ~State() = default;
    };

    /**
     * Destructor.
     */
    virtual ~SearchlightUDF() = default;

    /**
     * Save state to a buffer.
     *
     * @return pointer to the buffer or nullptr if nothing to save
     */
    virtual const State *SaveState() const = 0;

    /**
     * Load previously saved state.
     *
     * It is assumed the state was saved by the SaveState function.
     *
     * @param state saved state
     */
    virtual void LoadState(const State *state) = 0;

    /**
     * Returns the minimum of the expression among all windows.
     *
     * @return the minimum of the expression
     */
    virtual int64 Min() const override {
        if (!ComputeMinMax()) {
            adapter_->SetCustomFail();
            solver()->Fail();
            //return kint64max;
        }
        return udf_state_->min_;
    }

    /**
     * Sets the minimum of the expression.
     *
     * Currently, we do not propagate it to the variables, but we fail the
     * search if the required minimum is impossible.
     *
     * @param m the minimum to set
     */
    virtual void SetMin(int64 m) override {
        if (!ComputeMinMax()) {
            adapter_->SetCustomFail();
            solver()->Fail();
        } else if (m > udf_state_->max_) {
            solver()->Fail();
        }
    }

    /**
     * Returns the maximum of the expression among all windows.
     *
     * @return the maximum of the expression
     */
    virtual int64 Max() const override {
        if (!ComputeMinMax()) {
            adapter_->SetCustomFail();
            solver()->Fail();
            //return kint64min;
        }
        return udf_state_->max_;
    }

    /**
     * Sets the maximum of the expression.
     *
     * Currently, we do not propagate it to the variables, but we fail the
     * search if the required maximum is impossible.
     *
     * @param m the maximum to set
     */
    virtual void SetMax(int64 m) override {
        if (!ComputeMinMax()) {
            adapter_->SetCustomFail();
            solver()->Fail();
        } else if (m < udf_state_->min_) {
            solver()->Fail();
        }
    }

    /**
     * Attaches a demon to watch for range changes of the expression.
     *
     * @param d the demon to attach
     */
    virtual void WhenRange(Demon* d) override {
        for (auto var: orig_vars_) {
            var->WhenRange(d);
        }
    }

protected:
    /*
     * Constructor.
     *
     * @param s CP or-tools solver
     * @param state UDF state
     * @param adapter adapter to use for data access
     * @param attr parameter attribute
     * @param vars parameter variables
     */
    SearchlightUDF(Solver *s, State *state, const AdapterPtr &adapter,
            AttributeID attr, const std::vector<IntExpr *> &vars);

    /*
     * ParameterVar is a proxy to access variables and expressions.
     */
    class ParameterVar {
    public:
        // Constructor
        ParameterVar(const IntExpr *expr);

        // Is it bound?
        bool Bound() const {
            return min_ == max_;
        }

        // Value for bound expression
        int64 Value() const {
            return min_;
        }

        // Min
        int64 Min() const {
            return min_;
        }

        // Min
        int64 Max() const {
            return max_;
        }

        // Domain size
        uint64 Size() const;

        // Does it contain x?
        bool Contains(int64 x) const;

        // Debug string
        std::string DebugString() const {
            return expr_->DebugString();
        }

        // Return iterator.
        IntVarIterator *Iterator() const {
            return iter_;
        }

    private:
        friend class SearchlightUDF;

        // Range iterator for general IntExpr
        class RangeIterator : public IntVarIterator {
        public:
            // Constructor
            RangeIterator(const ParameterVar *var) : var_{var} {}

            // This method must be called before each loop.
            virtual void Init() override {
                min_ = var_->Min();
                max_ = var_->Max();
                v_ = min_;
                valid_ = min_ <= max_;
            }

            // This method indicates if we can call Value() or not.
            virtual bool Ok() const override {
                return valid_;
            }

            // This method returns the value of the hole.
            virtual int64 Value() const override {
                return v_;
            }

            // This method moves the iterator to the next value.
            virtual void Next() override {
                if (!valid_) return;
                if (v_ == max_) {
                    valid_ = false;
                } else {
                    ++v_;
                }
            }

            // Pretty Print.
            virtual std::string DebugString() const override {
                return "ParameterVar::RangeIterator";
            }

        private:
            // Var
            const ParameterVar *var_;
            // Valid
            bool valid_ = false;
            // Current val
            int64 v_ = 0, min_ = 0, max_ = 0;
        };

        // Update values
        void Update();

        // Is it real IntVar?
        const bool real_var_;
        // Min/max values
        int64 min_ = 0, max_ = 0;
        // Original expr
        const IntExpr *expr_;
        // Iterator
        IntVarIterator *iter_;
    };

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

    // Get current variables
    const std::vector<ParameterVar> &GetVars() const {
        // Might be new min/max; recache
        for (auto &v: param_vars_) {
            v.Update();
        }
        return param_vars_;
    }

    // Visit original IntExpr
    void VisitOriginalExprs(ModelVisitor *visitor) const;

    // Adapter
    const AdapterPtr adapter_; // the adapter for data access
    const AttributeID attr_; // parameter attribute

private:
    // Computes Min/Max
    virtual bool ComputeMinMax() const = 0;

    // Pointer to the UDF state (taken from the derived class)
    State * const udf_state_;
    // Original vars
    const std::vector<IntExpr *> orig_vars_;
    // Parameter variables
    mutable std::vector<ParameterVar> param_vars_;
};

/**
 * Saved state for UDF functions.
 */
using UDFStates = std::vector<std::unique_ptr<const SearchlightUDF::State>>;

} /* namespace searchlight */
#endif /* SEARCHLIGHT_SEARCHLIGHT_UDFS_H_ */

/* Copyright 2015, Brown University, Providence, RI.
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
 * @file relax.h
 *
 * This file defines stuff for the query relaxation subsystem.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_RELAX_H_
#define SEARCHLIGHT_RELAX_H_

#include "ortools_inc.h"
#include "searchlight.h"

namespace searchlight {

/**
 * Relaxable constraint. Represents a constraint that can be changed during the
 * search.
 */
class RelaxableConstraint : public Constraint {
public:
    /**
     * Relaxable constraint tags.
     */
    static const char BetweenConstTag[];
    static const char LessEqConstTag[];
    static const char GreaterEqConstTag[];

    /**
     * Construct a new relaxable constraint.
     *
     * @param solver solver instance
     */
    RelaxableConstraint(Solver* const solver) : Constraint(solver) {}
    /**
     * Sets new params for the constraint.
     *
     * The constraint is expected to be of interval type. For example,
     * arithmetic comparisons and between satify this property.
     *
     * @param l lower interval value
     * @param h high interval value
     */
    virtual void SetNewParams(int64 l, int64 h) = 0;

    /**
     * Get constraint params.
     *
     * @param l low interval value
     * @param h high interval value
     */
    virtual void GetParams(int64 &l, int64 &h) const = 0;

    /**
     * Return constraint's expression.
     *
     * @return constraint's expression
     */
    virtual IntExpr *GetExpr() const = 0;

    /**
     * Accept model visitor.
     *
     * @param visitor model visitor
     */
    virtual void Accept(ModelVisitor* const visitor) const override {
        visitor->VisitIntegerArgument(ModelVisitor::kValueArgument, id_);
    }

    /**
     * Set constraint id.
     *
     * @param id constraint id
     */
    void SetId(int64 id) {
        id_ = id;
    }

    /**
     * Return constraint ID.
     *
     * @return constraint ID
     */
    int64 Id() const {
        return id_;
    }

private:
	// Constraint id
	int64 id_ = -1;
};

/**
 * Relaxable between constraint.
 */
class BetweenCt : public RelaxableConstraint {
public:
    BetweenCt(Solver* const s, IntExpr* const v, int64 l, int64 u) :
        RelaxableConstraint(s),
        expr_(v),
        min_(l),
        max_(u) {}

    /**
     * Set new interval for the between constraint.
     *
     * @param l low interval value
     * @param h high interval value
     */
    virtual void SetNewParams(int64 l, int64 h) override {
        solver()->SaveAndSetValue(&min_, l);
        solver()->SaveAndSetValue(&max_, h);
    }

   /**
    * Get interval for the between constraint.
    *
    * @param l low interval value
    * @param h high interval value
    */
    virtual void GetParams(int64 &l, int64 &h) const override {
        l = min_;
        h = max_;
    }

    virtual void Post() {
        if (!expr_->IsVar()) {
            Demon* const d = solver()->
                    MakeConstraintInitialPropagateCallback(this);
            expr_->WhenRange(d);
        }
    }

    virtual void InitialPropagate() {
        expr_->SetRange(min_, max_);
    }

    virtual std::string DebugString() const {
        return StringPrintf(
                "RelaxableBetweenCt(%s, %" GG_LL_FORMAT "d, %" GG_LL_FORMAT "d)",
                expr_->DebugString().c_str(), min_, max_);
    }

    virtual void Accept(ModelVisitor* const visitor) const {
        visitor->BeginVisitConstraint(
                RelaxableConstraint::BetweenConstTag, this);
        RelaxableConstraint::Accept(visitor);
        visitor->VisitIntegerArgument(ModelVisitor::kMinArgument, min_);
        visitor->VisitIntegerExpressionArgument(
                ModelVisitor::kExpressionArgument, expr_);
        visitor->VisitIntegerArgument(ModelVisitor::kMaxArgument, max_);
        visitor->EndVisitConstraint(ModelVisitor::kBetween, this);
    }

    virtual IntExpr *GetExpr() const override {
        return expr_;
    }

private:
    IntExpr* const expr_;
    int64 min_;
    int64 max_;
};

/**
 * Relaxable >= constraint.
 */
class GreaterEqExprCst : public RelaxableConstraint {
public:
    GreaterEqExprCst(Solver* const s, IntExpr* const e, int64 v) :
        RelaxableConstraint(s),
        expr_(e),
        value_(v),
        demon_(nullptr) {}

    virtual void Post() {
        if (!expr_->IsVar() && expr_->Min() < value_) {
            demon_ = solver()->MakeConstraintInitialPropagateCallback(this);
            expr_->WhenRange(demon_);
        } else {
            // Clean the demon in case the constraint is posted during search.
            demon_ = nullptr;
        }
	}

    virtual IntExpr *GetExpr() const override {
        return expr_;
    }

    virtual void InitialPropagate() {
        expr_->SetMin(value_);
        if (demon_ != nullptr && expr_->Min() >= value_) {
            demon_->inhibit(solver());
        }
    }

    virtual std::string DebugString() const {
        return StringPrintf("(relaxable %s >= %" GG_LL_FORMAT "d)",
                            expr_->DebugString().c_str(), value_);
    }

    virtual IntVar* Var() {
        return solver()->MakeIsGreaterOrEqualCstVar(expr_->Var(), value_);
    }

    virtual void Accept(ModelVisitor* const visitor) const {
        visitor->BeginVisitConstraint(
                RelaxableConstraint::GreaterEqConstTag, this);
        RelaxableConstraint::Accept(visitor);
        visitor->VisitIntegerExpressionArgument(
                ModelVisitor::kExpressionArgument, expr_);
        visitor->VisitIntegerArgument(ModelVisitor::kValueArgument, value_);
        visitor->EndVisitConstraint(ModelVisitor::kGreaterOrEqual, this);
    }

    /**
     * Set new interval for the >= constraint.
     *
     * Note, h is actually ignore since it's always +inf for this constraint.
     *
     * @param l low interval value
     * @param h high interval value (ignored)
     */
    virtual void SetNewParams(int64 l, int64 h) override {
        solver()->SaveAndSetValue(&value_, l);
    }

    /**
     * Get interval for the >= constraint.
     *
     * High value is always +inf for this constraint.
     *
     * @param l low interval value
     * @param h high interval value
     */
    virtual void GetParams(int64 &l, int64 &h) const override {
	    l = value_;
	    h = kint64max;
    }

private:
    IntExpr* const expr_;
    int64 value_;
    Demon* demon_;
};

/**
 * Relaxable <= constraint.
 */
class LessEqExprCst : public RelaxableConstraint {
public:
    LessEqExprCst(Solver* const s, IntExpr* const e, int64 v) :
        RelaxableConstraint(s),
        expr_(e),
        value_(v),
        demon_(nullptr) {}

    virtual void Post() {
        if (!expr_->IsVar() && expr_->Max() > value_) {
            demon_ = solver()->MakeConstraintInitialPropagateCallback(this);
            expr_->WhenRange(demon_);
        } else {
            // Clean the demon in case the constraint is posted during search.
            demon_ = nullptr;
        }
    }

    virtual void InitialPropagate() {
        expr_->SetMax(value_);
        if (demon_ != nullptr && expr_->Max() <= value_) {
          demon_->inhibit(solver());
        }
    }

    virtual std::string DebugString() const {
        return StringPrintf("(relaxable %s <= %" GG_LL_FORMAT "d)",
                            expr_->DebugString().c_str(), value_);
    }

    virtual IntVar* Var() {
        return solver()->MakeIsLessOrEqualCstVar(expr_->Var(), value_);
    }

    virtual void Accept(ModelVisitor* const visitor) const {
        visitor->BeginVisitConstraint(
                RelaxableConstraint::LessEqConstTag, this);
        RelaxableConstraint::Accept(visitor);
        visitor->VisitIntegerExpressionArgument(
                ModelVisitor::kExpressionArgument, expr_);
        visitor->VisitIntegerArgument(ModelVisitor::kValueArgument, value_);
        visitor->EndVisitConstraint(ModelVisitor::kLessOrEqual, this);
    }

    virtual IntExpr *GetExpr() const override {
        return expr_;
    }

    /**
    * Set new interval for the <= constraint.
    *
    * Note, l is actually ignore since it's always -inf for this constraint.
    *
    * @param l low interval value (ignored)
    * @param h high interval value
    */
    virtual void SetNewParams(int64 l, int64 h) override {
        solver()->SaveAndSetValue(&value_, h);
    }

    /**
    * Get interval for the <= constraint.
    *
    * Low value is always -inf for this constraint.
    *
    * @param l low interval value
    * @param h high interval value
    */
    virtual void GetParams(int64 &l, int64 &h) const override {
        l = kint64min;
        h = value_;
    }

private:
    IntExpr* const expr_;
    int64 value_;
    Demon* demon_;
};

/**
 * This class represents constraint. Constraint is basically an interval,
 * describing all valid integer points.
 */
class IntervalConstraint {
public:
	/**
	 * Type of the constraint.
	 */
	enum class ConstraintType {
		LEQ,//!< LEQ Less than or equal
		GEQ,//!< GEQ Greater than or equal
		BET //!< BET Between
	};

	/**
	 * Create a new constraint.
	 *
	 * @param type constraint type
	 * @param min min value for the constraint interval
	 * @param max max value for the constraint interval
	 */
	IntervalConstraint(ConstraintType type, int64 min, int64 max) :
		type_(type),
		min_(min),
		max_(max) {

		if (type == ConstraintType::LEQ) {
			min_ = kint64min;
		} else if (type == ConstraintType::GEQ) {
			max_ = kint64max;
		}
	}

	/**
	 * Compute distance from the given point to the constraint interval.
	 *
	 * @param point point
	 * @return distance from point to constraint interval
	 */
	int64 Distance(int64 point) const {
		if (point >= min_ && point <= max_) {
			return 0;
		} else if (point > max_) {
			return point - max_;
		} else {
			return min_ - point;
		}
	}

	/**
	 * Compute the distance from this interval to the target interval.
	 *
	 * If the intervals intersect, the distance is supposed to be 0. Otherwise,
	 * we compute the minimum and maximum distance from the target interval's
	 * points to the nearest point of this interval.
	 *
	 * The return value indicates the relative position of the intervals.
	 * 1 means the target is on the right, -1 that it's on the left, 0 that
	 * they intersect.
	 *
	 * @param l low target interval point
	 * @param h high target interval point
	 * @param min_d minimum distance (returned)
	 * @param max_d maximum distance (returned)
	 * @return relative position of the intervals
	 */
	int MinMaxDistance(int64 l, int64 h, int64 &min_d, int64 &max_d) const {
		if (l > max_) {
			// Interval is to the right
			min_d = l - max_;
			max_d = h - max_;
			return 1;
		} else if (h < min_) {
			// Interval is to the left
			min_d = min_ - h;
			max_d = min_ - l;
			return -1;
		} else {
			// Intersect
			min_d = max_d = 0;
			return 0;
		}
	}

	/**
	 * Return constraint's minimum interval value.
	 *
	 * @return minimum interval value
	 */
	int64 Min() const {
		return min_;
	}

	/**
	 * Return constraint's maximum interval value.
	 *
	 * @return maximum interval value
	 */
	int64 Max() const {
		return max_;
	}

	/**
	 * Return constraint's type.
	 *
	 * @return constraint's type
	 */
	ConstraintType Type() const {
		return type_;
	}

private:
	// Constraint type
	ConstraintType type_;

	// Constraint interval
	int64 min_, max_;
};


/**
 * Relaxator is responsible for all the query relaxation logic.
 *
 * Each constraint is supposed to be some expression f(x) constrained over
 * an interval. It is also supposed to have a unique name across all solvers.
 *
 * The relaxation degree is: dw * dist + (1-dw)*viol_consts, where:
 * 	dw -- distance weight [0, 1] (query param)
 * 	dist -- relaxation distance [0, 1]
 * 	viol_consts -- violated constraints (relative, [0, 1])
 *
 */
class Relaxator : private boost::noncopyable {
public:
	/**
	 * Create a new relaxator instance.
	 *
	 * @param sl Searchlight instance
	 * @param solvers
	 * @param dist_w distance weight for the relaxation degree
	 * @param res_num the number of results to track (top-k)
	 */
	Relaxator(Searchlight &sl, size_t solvers, double dist_w, size_t res_num) :
		sl_(sl),
		solvers_num_(solvers),
		distance_weight_(dist_w),
		res_num_(res_num) {}

	/**
	 * Register new relaxation constraint for a particular solver.
	 *
	 * @param name constraint name
	 * @param solver_id solver id
	 * @param constr constraint object
	 * @param expr constraint's expression object
	 * @param max_l maximum relaxation for l
	 * @param max_h maximum relaxation for h
	 */
	void RegisterConstraint(const std::string &name, size_t solver_id,
			RelaxableConstraint *constr, IntExpr *expr,
			int64 max_l, int64 max_h);

	/**
	 * Register a solver fail from a particular Searchlight solver.
	 *
	 * Check if the fail can be relaxed in the future, and if so, then save it.
	 * Otherwise, just ignore it.
	 *
	 * @param solver_id SearchlightSolver id
	 */
	void RegisterFail(size_t solver_id);

	/**
	 * Update Lower Relaxation Degree (LDR) bound.
	 *
	 * @param lrd new lower relaxation degree bound
	 */
	void UpdateLRD(double lrd) {
		assert(lrd >= 0.0 && lrd <= 1.0);
		// Relaxed is fine here; we don't need any specific ordering
		lrd_.store(lrd, std::memory_order_relaxed);
	}

	/**
	 * Return current LRD.
	 *
	 * @return current LRD
	 */
	double GetLRD() const {
		return lrd_.load(std::memory_order_relaxed);
	}

	/**
	 * Check if the query is a relaxing one.
	 *
	 * @return true, if we're relaxing the query; false otherwise
	 */
	bool RelaxingQuery() const {
		return !orig_consts_.empty();
	}

	void AddResult(const std::vector<int64> &const_expr_vals);

private:
	/*
	 * Compute the maximum [0, 1] relaxation distance for the specified number
	 * of violated constraints.
	 */
	double MaxUnitRelaxDistance(size_t viol_constrs) const {
		return lrd_.load(std::memory_order_relaxed) -
				(1 - distance_weight_) * double(viol_constrs) /
				orig_consts_.size() / distance_weight_;
	}

	// Compute relaxation distance
	double RelaxDistance(double relax_dist, size_t viol_const_num) const {
		return distance_weight_ * relax_dist +
				(1 - distance_weight_) *
				double(viol_const_num) / orig_consts_.size();
	}

	// Information about each original constraint
	struct ConstraintInfo {
		// Original constraint interval
		const IntervalConstraint int_;

		// Maximum relaxation
		const int64 max_l_, max_h_;

		// Pointers to each constraint object for each solver
		std::vector<RelaxableConstraint *> solver_const_;

		// Pointers to each constraint's expression object for each solver
		std::vector<const IntExpr *> solver_exprs_;

		// Construct a new constraint info.
		ConstraintInfo(const IntervalConstraint &constr, size_t solver_num,
				int64 max_l, int64 max_h) :
			int_{constr},
			max_l_{max_l},
			max_h_{max_h},
			solver_const_(solver_num),
			solver_exprs_(solver_num) {}

		// Return maximum relaxation distance to the left or right
		int64 MaxRelaxDist(bool left) const {
			if (left) {
				if (max_l_ == kint64min) {
					return kint64min;
				} else {
					return int_.Min() - max_l_;
				}
			} else {
				if (max_h_ == kint64max) {
					return kint64max;
				} else {
					return max_h_ - int_.Max();
				}
			}
		}
	};

	// Saved fail for later replay
	struct FailReplay {
		struct FailedConstraint {
			// Constraint id
			size_t const_id_;

			// Low/high values at the fail point
			int64 l_, h_;

			// Computed relaxation interval
			int64 rl_, rh_;
		};

		// Failed constraints
		std::vector<FailedConstraint> failed_const_;

		// Saved decision variables for the replay
		std::vector<IntVarDomainInfo> saved_vars_;

		// Relaxation degree
		double relax_degree_;

		// relax degree-based comparison
		bool operator>(const FailReplay &other) const {
			return relax_degree_ > other.relax_degree_;
		}
	};

	// Queue of replays ranked by relaxation degree (lowest first)
	std::priority_queue<FailReplay, std::vector<FailReplay>,
		std::greater<FailReplay>> fail_replays_;

	// Fill in new relaxation interval; return relaxation distance (-1 failed)
	int64 ComputeNewReplayInterval(FailReplay::FailedConstraint &failed_const,
			double max_relax_dist, int rel_pos);

	// Main Searchlight reference
	Searchlight &sl_;

	// Total number of solvers
	const size_t solvers_num_;

	// Mapping from the constraint's name to the constraint's id.
	std::unordered_map<std::string, size_t> const_to_id_;

	// Original constraints
	std::vector<ConstraintInfo> orig_consts_;

	// Lower relaxation degree, [0, 1]
	std::atomic<double> lrd_{1.0};

	// Distance weight for the relaxation degree (query param)
	const double distance_weight_;

	// Maximum number of results to track
	const size_t res_num_;

	// Map of top results by relaxation degree
	std::map<double, size_t> top_res_rd_;

	// For concurrency control
	std::mutex mtx_;
};

/**
 * This search monitor catches fails and calls relaxator to record them for
 * possible future replays.
 */
class FailCollectorMonitor : public SearchMonitor {
public:
	/**
	 * Create a new FailCollectorMonitor.
	 *
	 * @param solver CP solver
	 * @param rel query relaxator
	 * @param sid Searchlight solver id
	 */
	FailCollectorMonitor(Solver &solver, Relaxator &rel, size_t sid) :
		SearchMonitor(&solver),
		relaxator_(rel),
		solver_id_(sid) {}

    /**
     * Called at the beginning of a fail.
     *
     * Just call the relaxator. It knows what to do.
     */
	virtual void BeginFail() override {
    	relaxator_.RegisterFail(solver_id_);
    }

private:
	// Query relaxator
	Relaxator &relaxator_;

	// Searchlight solver id
	size_t solver_id_;
};
} /* namespace searchlight */
#endif /* SEARCHLIGHT_RELAX_H_ */

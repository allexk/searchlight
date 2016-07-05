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
#include "searchlight_udfs.h"

#include <queue>
#include <mutex>
#include <atomic>

namespace searchlight {

class Searchlight;

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
     * Model parameter tags for Protobuf import/export.
     */
    static const char ModelIDTag[];

    /**
     * Check if the constraint is relaxable.
     *
     * The check is based on the type name.
     *
     * @param type type name
     * @return true, if the constraint is relaxable; false, otherwise
     */
    static bool IsRelaxable(const std::string &type) {
        return type == BetweenConstTag ||
               type == LessEqConstTag ||
               type == GreaterEqConstTag;
    }

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
        visitor->VisitIntegerArgument(ModelIDTag, id_);
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
 * Relaxable constraints found in the model.
 *
 * Ordered by IDs, which are issued by the relaxator before serialization.
 */
using RelaxableConstraints = std::vector<RelaxableConstraint *>;

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

    virtual void Post() override {
        if (!expr_->IsVar()) {
            Demon* const d = solver()->
                    MakeConstraintInitialPropagateCallback(this);
            expr_->WhenRange(d);
        }
    }

    virtual void InitialPropagate() override {
        expr_->SetRange(min_, max_);
    }

    virtual std::string DebugString() const override {
        return StringPrintf(
                "RelaxableBetweenCt(%s, %" GG_LL_FORMAT "d, %" GG_LL_FORMAT "d)",
                expr_->DebugString().c_str(), min_, max_);
    }

    virtual void Accept(ModelVisitor* const visitor) const override {
        visitor->BeginVisitConstraint(
                RelaxableConstraint::BetweenConstTag, this);
        RelaxableConstraint::Accept(visitor);
        visitor->VisitIntegerArgument(ModelVisitor::kMinArgument, min_);
        visitor->VisitIntegerExpressionArgument(
                ModelVisitor::kExpressionArgument, expr_);
        visitor->VisitIntegerArgument(ModelVisitor::kMaxArgument, max_);
        visitor->EndVisitConstraint(RelaxableConstraint::BetweenConstTag, this);
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

    virtual void Post() override {
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

    virtual void InitialPropagate() override {
        expr_->SetMin(value_);
        if (demon_ != nullptr && expr_->Min() >= value_) {
            demon_->inhibit(solver());
        }
    }

    virtual std::string DebugString() const override {
        return StringPrintf("(relaxable %s >= %" GG_LL_FORMAT "d)",
                            expr_->DebugString().c_str(), value_);
    }

    virtual void Accept(ModelVisitor* const visitor) const override {
        visitor->BeginVisitConstraint(
                RelaxableConstraint::GreaterEqConstTag, this);
        RelaxableConstraint::Accept(visitor);
        visitor->VisitIntegerExpressionArgument(
                ModelVisitor::kExpressionArgument, expr_);
        visitor->VisitIntegerArgument(ModelVisitor::kValueArgument, value_);
        visitor->EndVisitConstraint(
                RelaxableConstraint::GreaterEqConstTag, this);
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

    virtual void Post() override {
        if (!expr_->IsVar() && expr_->Max() > value_) {
            demon_ = solver()->MakeConstraintInitialPropagateCallback(this);
            expr_->WhenRange(demon_);
        } else {
            // Clean the demon in case the constraint is posted during search.
            demon_ = nullptr;
        }
    }

    virtual void InitialPropagate() override {
        expr_->SetMax(value_);
        if (demon_ != nullptr && expr_->Max() <= value_) {
          demon_->inhibit(solver());
        }
    }

    virtual std::string DebugString() const override {
        return StringPrintf("(relaxable %s <= %" GG_LL_FORMAT "d)",
                            expr_->DebugString().c_str(), value_);
    }

    virtual void Accept(ModelVisitor* const visitor) const override {
        visitor->BeginVisitConstraint(
                RelaxableConstraint::LessEqConstTag, this);
        RelaxableConstraint::Accept(visitor);
        visitor->VisitIntegerExpressionArgument(
                ModelVisitor::kExpressionArgument, expr_);
        visitor->VisitIntegerArgument(ModelVisitor::kValueArgument, value_);
        visitor->EndVisitConstraint(RelaxableConstraint::LessEqConstTag, this);
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
	 * Return relative position of the interval to this one.
	 *
	 * @param l left interval bound
	 * @param h right interval bound
	 * @return -1 -- to the left; 1 -- to the right; 0 -- intersect
	 */
	int RelIntervalPos(int64 l, int64 h) const {
	    return l > max_? 1 : (h < min_ ? -1 : 0);
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
 * Results contractor.
 */
class Contractor {
public:
    /**
     * Ctor.
     *
     * @param relaxator relaxator
     * @param card cardinality
     * @param constrs constraints to contract
     */
    Contractor(Relaxator &relaxator, size_t card,
               const std::vector<size_t> &constrs) :
                   relaxator_(relaxator),
                   card_req_(card),
                   contr_constrs_(constrs) {}
    /**
     * Dtor.
     */
    virtual ~Contractor() {}

    /**
     * Check if the solution is valid.
     *
     * @param sol solution
     * @param update true, if we can update the results; false, just check
     * @return true, if valid
     */
    virtual bool CheckSolution(const LiteVarAssignment &sol, bool update) = 0;

    /**
     * Register new solution.
     *
     * @param sol solution
     * @param rank rank
     */
    virtual void RegisterSolution(const std::vector<int64_t> &sol,
                                  double rank) = 0;

protected:
    // Relaxator
    Relaxator &relaxator_;
    // Cardinality reqs
    size_t card_req_;
    // Constraints we contract
    std::vector<size_t> contr_constrs_;
    // Mutex
    std::mutex mtx_;
};

class SkylineContractor : public Contractor {
public:
    SkylineContractor(Relaxator &relaxator, size_t card,
               const std::vector<size_t> &constrs) :
                   Contractor(relaxator, card, constrs) {}

    virtual bool CheckSolution(const LiteVarAssignment &sol,
                               bool update) override {

    }

    virtual void RegisterSolution(const std::vector<int64_t> &sol,
                                  double rank) override {

    }
private:
    // Dominating solutions
    std::list<std::vector<int64_t>> sols_;
};

class RankContractor : public Contractor {
public:
    RankContractor(Relaxator &relaxator, size_t card,
               const std::vector<size_t> &constrs,
               const std::vector<bool> &spec) :
                   Contractor(relaxator, card, constrs),
                   spec_(spec) {
        // Checking
        assert(constrs.size() == spec.size());
    }

    virtual bool CheckSolution(const LiteVarAssignment &sol,
                               bool update) const override {

    }

    virtual void RegisterSolution(const std::vector<int64_t> &sol,
                                  double rank) override {

    }

private:
    // true/false maximize for each constraint
    std::vector<bool> spec_;
    // Rank priority queue
    std::priority_queue<double, std::vector<double>, std::greater<double>> ranks_;
    // Current top lower bound rank
    std::atomic<double> lr_{0.0};
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
     * Relaxator heuristic when registering a fail.
     */
    enum class RegisterHeuristic {
        GUESS, ///!< GUESS Tries to guess the violated constraint.
        GUESS_ALL, ///! GUESS_ALL Guess for initial, ALL for repeat fails.
        ALL   ///!< ALL Checks all the constraints for violations.
    };

    /**
     * Relaxation type to apply to replays.
     */
    enum class ReplayRelaxation {
        VIOLATED, ///!< VIOLATED Relax only violated constraints.
        ALL       ///!< ALL Relax all constraints.
    };

public:
	/**
	 * Create a new relaxator instance.
	 *
	 * @param sl Searchlight instance
	 * @param solvers number of solvers
	 * @param dist_w distance weight for the relaxation degree
	 * @param res_num the number of results to track (top-k)
	 */
	Relaxator(Searchlight &sl, size_t solvers, double dist_w, size_t res_num);

	/**
	 * Destructor.
	 */
	~Relaxator();

	/**
	 * Register new relaxation constraint for a particular solver.
	 *
	 * @param name constraint name
	 * @param solver_id solver id
	 * @param constr constraint object
	 * @param max_l maximum relaxation for low bound
	 * @param max_h maximum relaxation for high bound
	 */
	void RegisterConstraint(const std::string &name, size_t solver_id,
			RelaxableConstraint *constr, int64 max_l, int64 max_h);

	/**
	 * Register a solver fail from a particular Searchlight solver.
	 *
	 * Check if the fail can be relaxed in the future, and if so, then save it.
	 * Otherwise, just ignore it.
	 *
	 * @param solver_id SearchlightSolver id
     * @param imm_fail true, if the fail is immediate (no new decisions)
	 */
	void RegisterFail(size_t solver_id, bool imm_fail) {
	    RegisterHeuristic h = register_heur_;
	    if (h == RegisterHeuristic::GUESS_ALL) {
	        h = solver_info_[solver_id].in_replay_ ? RegisterHeuristic::ALL :
	                RegisterHeuristic::GUESS;
	    }
	    RegisterFail(solver_id, h, imm_fail);
	}

	/**
	 * Compute current VC-specification and best RD for it.
	 *
	 * Current means for the current state of the solver. The function should
	 * be called when the solver is paused, e.g., at a leaf.
	 *
	 * @param sid local solver id
	 * @param rd best relaxation degree (out)
	 * @param vc VC-specification (out)
	 * @return true, if the current state is relaxable; fals, otherwise
	 */
	bool ComputeCurrentVCAndRD(size_t sid, double &rd, Int64Vector &vc) const;

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

	/**
	 * Compute the relaxation degree of the result.
	 *
	 * The result parameter is supposed to contain resulting values for all
	 * relaxable constraints. This function checks if the constraint is really
	 * violated.
	 *
	 * @param rel_result result for relaxable constraints
	 * @return result's relaxation degree
	 */
	double ComputeResultRelaxationDegree(const Int64Vector &rel_result) const;

	/**
	 * Compute the best relaxation degree of the vc specification.
	 *
	 * Since vc specification contains relaxation intervals for the constraints,
	 * we can compute the best and worst case relaxations. This function
	 * returns the best one, based on the provided intervals.
	 *
	 * @param vc violated constraint specification
	 * @return relaxation degree for vc
	 */
	double ComputeViolSpecBestRelaxationDegree(const Int64Vector &vc) const;

	/**
	 * Report a new result with the specified relaxation degree.
	 *
	 * @param rd relaxation degree
	 */
	void ReportResult(double rd);

	/**
	 * Report a new solution/rank.
	 *
	 * @param sol solution
	 * @param rank rank
	 */
	void ReportRankSol(const std::vector<int64_t> &sol, double rank) const {
	    contractor_->RegisterSolution(sol, rank);
	}

	/**
	 * Check if there are any fail replays left
	 *
	 * @return true, if replays left; false, otherwise
	 */
	bool HasReplays() const {
	    std::lock_guard<std::mutex> lock{mtx_};
	    return !fail_replays_.empty() || !forced_replays_queue_.empty();
	}

	/**
	 * Apply violated constraints specification to the solver.
	 *
	 * @param vc_spec violated constraint specification
	 * @param solver_id solver id
	 */
	void ApplyViolatedConstSpec(const Int64Vector &vc_spec,
	                            size_t solver_id) const {
	    for (size_t i = 0; i < vc_spec.size(); i += 3) {
	        const size_t const_id = vc_spec[i];
	        const int64_t l = vc_spec[i + 1];
	        const int64_t h = vc_spec[i + 2];
	        orig_consts_[const_id].solver_const_[solver_id]->SetNewParams(l, h);
	    }
	}

	/**
	 * Return the next fail replay, if available.
	 *
	 * @param solver_id local solver id
	 * @param dom_info variable domains info
	 * @param vc_spec violated constraints spec
	 * @param udfs saved udf state (out)
	 * @return true, if there was a replay; false, otherwise
	 */
	bool GetFailReplay(size_t solver_id,
	        std::vector<IntVarDomainInfo> &dom_info, Int64Vector &vc_spec,
	        UDFStates &udfs);

	/**
	 * Compute VC specification with the maximum relaxation possible.
	 *
	 * The maximum relaxation assumes all constraints are violated and
	 * corresponds to the maximum possible relaxation distances for each
	 * according to the LRD.
	 *
	 * @param viol_num number of violated constraints
	 * @return VC specification with maximal relaxation
	 */
	Int64Vector GetMaximumRelaxationVCSpec(size_t viol_num) const;

	/**
	 * Check if re-replays might be needed.
	 *
	 * Re-replays happen when a replay fails again. With some replay methods,
	 * like ALL, this is not needed, since they already give the maximum
	 * possible relaxation.
	 *
	 * @return true, if re-replays might be needed; false, otherwise
	 */
	bool ReReplaysNeeded() const {
	    return replay_relax_ != ReplayRelaxation::ALL || replay_rd_ < 1.0;
	}

	/**
	 * Callback to call when the solver finished a replay.
	 *
	 * @param solver_id solver local id
	 */
	void ReplayFinished(size_t solver_id) {
	    assert(solver_info_[solver_id].in_replay_);
	    solver_info_[solver_id].in_replay_ = false;
	    solver_info_[solver_id].replay_.failed_consts_.clear();
	}

private:
	/**
	 * Replay sorting method.
	 *
	 * The method determines which replays will be replayed first later.
	 */
	enum class ReplaySortMethod {
	    BEST,   ///!< BEST sort based on best relaxation degree
	    WORST,  ///!< WORST sort based on worst relaxation degree
	    PROD,   ///!< PROD sort based on the product of relaxation degrees
	    TIME    ///!< TIME sort based on time registered.
	};

	/*
	 * Compute the maximum [0, 1] relaxation distance for the specified number
	 * of violated constraints.
	 */
	double MaxUnitRelaxDistance(size_t viol_constrs) const {
	    assert(viol_constrs);
	    const double res = (lrd_.load(std::memory_order_relaxed) -
				(1 - distance_weight_) * double(viol_constrs) /
				orig_consts_.size()) / distance_weight_;
	    /*
	     * The value might be greater than 1.0 in case the number of violated
	     * constraints is small and the LRD is still large. Still, the
	     * distance should not be greater than 1.0 -- the maximum one.
	     */
	    return std::min(res, 1.0);
	}

	// Compute relaxation distance
	double RelaxDistance(double relax_dist, size_t viol_const_num) const {
	    assert(relax_dist >= 0.0 && relax_dist <= 1.0);
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

		// Normalization interval for this function values
		const int64 min_norm_, max_norm_, max_norm_dist_;

        // Pointers to solver constraints
		std::vector<RelaxableConstraint *> solver_const_;

		// Construct a new constraint info.
		ConstraintInfo(const IntervalConstraint &constr, size_t solver_num,
				int64 max_l, int64 max_h) :
			int_{constr},
			max_l_{max_l},
			max_h_{max_h},
			min_norm_{constr.Min() == kint64min ? max_l : constr.Min()},
            max_norm_{constr.Max() == kint64max ? max_h : constr.Max()},
            max_norm_dist_{max_norm_ - min_norm_},
			solver_const_(solver_num) {}

		// True, if we can relax up to dist to the required direction
        bool CanRelax(bool left, int64_t dist) const {
            return MaxRelaxDist(left) >= dist;
		}

		// Return maximum relaxation distance to the left or right
		int64 MaxRelaxDist(bool left) const {
			if (left) {
				if (int_.Min() == kint64min || max_l_ == kint64min) {
					return kint64max;
				} else {
					return int_.Min() - max_l_;
				}
			} else {
				if (int_.Max() == kint64max || max_h_ == kint64max) {
					return kint64max;
				} else {
					return max_h_ - int_.Max();
				}
			}
		}

		// Compute maximum relaxation interval with respect to max_relax degree
		void MaxRelaxDist(double max_relax, int64_t &l, int64_t &h) const {
		    assert(max_relax >= 0.0 && max_relax <= 1.0);
		    // Left point
		    l = MaxRelaxDist(true);
		    if (l != kint64max) {
		        l = int_.Min() - l * max_relax;
		    }
		    // Right point
		    h = MaxRelaxDist(false);
		    if (h != kint64max) {
                h = int_.Max() + h * max_relax;
		    }
		}

		// Corrects maximum relaxation for the interval
		void CorrectMaxRelaxInerval(double max_relax, int64_t &l,
		                            int64_t &h) const {
		    int64_t min_l, max_h;
		    MaxRelaxDist(max_relax, min_l, max_h);
		    if (l > int_.Max()) {
		        // l-h is to the right
		        h = std::min(h, max_h);
		    } else if( h < int_.Min()) {
		        // l-h is to the left
		        l = std::max(l, min_l);
		    }
		}

		// Compute relative relaxation distance based on the point
		double RelRelaxDist(int64 x) const {
		    const int64 imin = int_.Min();
            const int64 imax = int_.Max();
		    if (x > imax) {
		        return double(x - imax) / MaxRelaxDist(false);
		    } else if (x < imin) {
		        return double(imin - x) / MaxRelaxDist(true);
		    } else {
		        return 0;
		    }
		}

		// Compute minimum relative relaxation distance for the interval
		double MinRelRelaxDist(int64 l, int64 h) const {
            const int64 imin = int_.Min();
            const int64 imax = int_.Max();
            /*
             * We assume that the interval is either completely to the left or
             * right of the original interval, since that is how we generate
             * them.
             */
            if (l > imax) {
                return double(l - imax) / MaxRelaxDist(false);
            } else if (h < imin) {
                return double(imin - h) / MaxRelaxDist(true);
            } else {
                return 0;
            }
		}

		// Compute rank of the point.
		double PointRank(int64 p, bool maximize) const {
		    if (p > max_norm_) {
		        p = max_norm_;
		    } else if (p < min_norm_) {
		        p = min_norm_;
		    }
		    double res = double(max_norm_ - p) / max_norm_dist_;
		    assert(res >= 0.0 && res <= 1.0);
		    if (maximize) {
		        res = 1.0 - res;
		    }
		    return res;
		}

		// Compute interval ranks
		void IntervalRank(int64 l, int64 h, bool maximize, double &minr,
		                   double &maxr) const {
		    minr = PointRank(l, maximize);
            maxr = PointRank(h, maximize);
            if (!maximize) {
                std::swap(minr, maxr);
            }
		}
	};

	// Saved fail for later replay
	struct FailReplay {
	    // Pointer
	    using Ptr = std::unique_ptr<FailReplay>;
	    // Info about a failed constraint
		struct FailedConstraint {
			// Constraint id
			size_t const_id_;

			// Left/right relaxation
			int rel_pos_;

			// Computed relaxation distance (low/high)
			int64 rl_, rh_;
		};

		// Failed constraints
		std::vector<FailedConstraint> failed_const_;

		// Saved decision variables for the replay
		std::vector<IntVarDomainInfo> saved_vars_;

		// Saved UDF states
		UDFStates saved_udfs_;

		// Relaxation degree: best and the worst one possible for the replay
		double best_relax_degree_, worst_relax_degree_;

		// Fail timestamp
		size_t timestamp_;
	};
	// To output the structure
	friend std::ostream &operator<<(std::ostream &os,
	    const Relaxator::FailReplay &fr);

    // Replay sorter
	struct ReplaySort {
	    // Constructs sorter based on method
	    ReplaySort(ReplaySortMethod method = ReplaySortMethod::BEST) :
	        sort_method_(method) {}

	    // Call operator for sorting
        bool operator()(const FailReplay &f1, const FailReplay &f2) const {
            switch (sort_method_) {
                case ReplaySortMethod::BEST:
                    return f1.best_relax_degree_ > f2.best_relax_degree_;
                case ReplaySortMethod::WORST:
                    return f1.worst_relax_degree_ > f2.worst_relax_degree_;
                case ReplaySortMethod::PROD:
                    return f1.best_relax_degree_ * f1.worst_relax_degree_ >
                        f2.best_relax_degree_ * f2.worst_relax_degree_;
                case ReplaySortMethod::TIME:
                    return f1.timestamp_ > f2.timestamp_;
                default:
                    assert(false);
                    return true;
            }
        }
	private:
	    ReplaySortMethod sort_method_;
    };

	/*
	 * Various stats for the Relaxator.
	 */
	struct RelaxatorStats {
	    // Total time spent for fail catching
	    std::chrono::microseconds total_fail_time_{0};
        // Total time spent for successful fail catching
        std::chrono::microseconds total_success_fail_time_{0};
	    // Total fails caught
	    std::atomic<size_t> total_fails_caught_{0};
        // Total fails registered
	    size_t total_fails_registered_{0};
        // Total fails replayed
	    size_t total_fails_replayed_{0};
	    // Total constraints taken from last replay on GUESS
	    std::atomic<size_t> total_const_reguessed_{0};
	    // Total fails registration retried with ALL heuristic
	    std::atomic<size_t> total_fails_heur_retried_{0};
	};
	// Output to stream
	friend std::ostream &operator<<(std::ostream &os,
	    const Relaxator::RelaxatorStats &rs);

    /*
     * Relaxator solver info.
     */
    struct SolverReplayInfo {
        // Last replay
        struct Replay {
            // Failed constraints: const_id -> FailedConstraint
            std::unordered_map<size_t, FailReplay::FailedConstraint> failed_consts_;
        };
        Replay replay_;

        // True, if it's replaying
        bool in_replay_ = false;

        // Frame to write stats to
        size_t GetStatsFrame() const {
            return in_replay_ ? 1 : 0;
        }
    };

    /*
     * Return violated constraint spec vector.
     *
     * The vector has the following format:
     *  <constraint id, left bound, right bound> triplet for each constraint
     */
	Int64Vector ViolConstSpec(const FailReplay &replay) const;

	// Internal register fail function
	void RegisterFail(size_t solver_id, RegisterHeuristic rh, bool imm_fail);

	// Fill in new relaxation interval; return true if possible to relax
	bool ComputeNewReplayInterval(FailReplay &replay,
	        FailReplay::FailedConstraint &failed_const,
			double max_relax_dist) const;

	// Compute relaxation for replay (true, if possible to relax)
	bool ComputeFailReplayRelaxation(FailReplay &replay) const;

	// Update time stats using the provided start time
	void UpdateTimeStats(RelaxatorStats &stats,
	    const std::chrono::steady_clock::time_point &start_time,
	    bool success);

	// Main Searchlight reference
	Searchlight &sl_;

    // Queue of replays ranked by relaxation degree (lowest first)
    std::priority_queue<FailReplay,
        std::vector<FailReplay>, ReplaySort> fail_replays_;

    // Forced fail replays (have priotorty over fail_replays_)
    std::queue<FailReplay> forced_replays_queue_;

	// Total number of solvers
	const size_t solvers_num_;

	// Info for each solver
	std::vector<SolverReplayInfo> solver_info_;

	// Mapping from the constraint's name to the constraint's id.
	std::unordered_map<std::string, size_t> const_to_id_;

	// Original constraints
	std::vector<ConstraintInfo> orig_consts_;

	// Lower relaxation degree, [0, 1]
	std::atomic<double> lrd_{1.0};

	// Fail timestamp
	size_t fail_stamp_{0};

	// Distance weight for the relaxation degree (query param)
	const double distance_weight_;

	// Relaxation degree to use with replays
	double replay_rd_;

	// Maximum number of results to track
	const size_t res_num_;

	// Do we save UDFs for replays?
	bool save_udfs_for_replay_;

	// Do we force some replays before others?
	bool force_replays_;

    // Max heap to count result relaxation degrees
    std::priority_queue<double> top_results_;

    // Relaxator stats: for several pahases
    RelaxatorStats stats_[2];

    // Register fail heuristic
    RegisterHeuristic register_heur_;

    // Replay relaxation type
    ReplayRelaxation replay_relax_;

    // Replay sort method
    ReplaySortMethod sort_method_;

    // Contractor
    std::unique_ptr<Contractor> contractor_;

	// For concurrency control
	mutable std::mutex mtx_;
};

class SearchlightSolver;


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
	 * @param sl_solver Searchlight solver
	 */
	FailCollectorMonitor(Solver &solver, const SearchlightSolver &sl_solver,
	        Relaxator &rel);

    /**
     * Called at the beginning of a fail.
     *
     * Just call the relaxator. It knows what to do.
     */
	virtual void BeginFail() override;

    /**
     * Called before the initial propagation.
     */
    virtual void BeginInitialPropagation() override {
        init_propagation_finished_ = false;
    }

    /**
     *  Called after the initial propagation.
     */
    virtual void EndInitialPropagation() override {
        init_propagation_finished_ = true;
    }

private:
	// Query relaxator
	Relaxator &relaxator_;

	// Searchlight solver instance
	const SearchlightSolver &sl_solver_;

	// Searchlight solver id
	const size_t solver_id_;

	// True, if the search went beyond initial propagation
	bool init_propagation_finished_ = false;
};
} /* namespace searchlight */
#endif /* SEARCHLIGHT_RELAX_H_ */

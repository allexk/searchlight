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
 * @file validator.cpp
 * The implementation of the validator.
 *
 * @author Alexander Kalinin
 */

#include "validator.h"
#include "base/callback.h"

namespace searchlight {

/**
 * A decision that just restores the given assignment causing the corresponding
 * variables to take the values from the assignment. In its right branch it
 * soes nothing, thus creating a right-deep tree.
 */
class RestoreAssignment : public Decision {
public:
    /**
     * Constructs a restore assignment decision.
     *
     * @param asgn the assignment to restore
     * @param callback_flag_ref pointer to a flag that should be set when
     *  Apply() is called
     */
    RestoreAssignment(AssignmentPtr asgn, bool *callback_flag_ref) :
        asgn_(asgn),
        callback_flag_ref_(callback_flag_ref) {}

    /**
     * The destructor.
     */
    virtual ~RestoreAssignment() {}

    /**
     * The left branch. In this case it restores the assignment.
     *
     * @param s the solver
     */
    virtual void Apply(Solver* const s) {
        // callback_flag should be set back after we backtrack
        s->SaveAndSetValue(callback_flag_ref_, true);
        asgn_->Restore();
    }

    /**
     * The right branch. It does nothing.
     *
     * @param s the solver
     */
    virtual void Refute(Solver* const s) {}

    /**
     * Returns a string for debug logging.
     *
     * @return a debug string
     */
    virtual std::string DebugString() const {
        return "RestoreAssignmentDecision";
    }

private:
    // The assignment to restore
    AssignmentPtr asgn_;

    // The flag to set when the Apply is calles
    bool *callback_flag_ref_;
};

/**
 * This is a DB that just generates decisions to restore the current assignment,
 * causing the corresponding variables to take the assignment's values.
 *
 * The main idea is that the validator is going to run this DB in a thread
 * until it is stopped. Then it will produce NULL, which will fail the search.
 * For this, we generate a right-deep tree: left node checks an assignment,
 * the right one does nothing.
 *
 * Since we have to produce a leaf after every successful restoration, this
 * db produces alternate decisions: after a successful restore a NULL is
 * generated that cuts the left branch. This is mostly taken care of by
 * the solver itself.
 *
 */
class RestoreAssignmentBuilder : public DecisionBuilder {
public:

    /**
     * Creates a new restore assignment builder.
     *
     * @param validator the validator producing assignments
     */
    explicit RestoreAssignmentBuilder(Validator &validator)
        : validator_(validator),
          just_restored_(false) {}

    /**
     * Destructor.
     */
    virtual ~RestoreAssignmentBuilder() {}


    /**
     * Sets variable values according to the assignment. Always returns
     * NULL, since search does not go beyond assigning values.
     *
     * @param solver the current solver
     * @return NULL
     */
    virtual Decision* Next(Solver* const solver) {
        if (just_restored_) {
           /*
            * We are here since the restoration was successful: cut the branch.
            * just_restored_ will fall back to 'false' automatically.
            */
            return NULL;
        }

        // check for termination
        if (validator_.CheckTerminate()) {
            // this will stop the search by cutting the right branch
            return NULL;
        }

        // pick the next assignment
        if (asgns_.empty()) {
            // need to ask the validator
            AssignmentPtrVector *new_asgns = validator_.GetNextAssignments();
            if (!new_asgns) {
                // No more assignments: stop the validator search
                return NULL;
            }
            asgns_.insert(asgns_.end(), new_asgns->begin(), new_asgns->end());
            delete new_asgns;
        }

        AssignmentPtr *next_asgn = asgns_.back();
        asgns_.pop_back();
        return solver->RevAlloc(new RestoreAssignment(next_asgn));
    }

    /**
     * Returns the name for debugging logging.
     * @return
     */
    virtual string DebugString() const {
        return "RestoreAssignment";
    }

 private:
    // The validator producing the assignments
    Validator validator_;

    // Assignments to validate
    AssignmentPtrVector asgns_;

    /*
     * Was the last decision a restoration? This flag is set by the decision
     * and is automatically restored during backtracking(via a fail and
     * leaf-based).
     */
    bool just_restored_;
};

Validator::Validator(const Searchlight &sl, const StringVector &var_names,
        SolutionCollector &collector) :
        sl_(sl),
        solver_("validator solver"),
        adapter_(sl.CreateAdapter()),
        search_vars_prototype_(&solver_),
        collector_(collector),
        terminate_(false),
        solver_status_(false) {
    /*
     * There might be a better way to do this, but... Since the solver does
     * not have a copying constructor, we first export the model into a
     * Protobuf buffer and then load it again via the loader, thus
     * duplicating the original solver.
     *
     * Note, no original monitors or decision builders are copied since
     * the search has not started yet.
     */
    CPModelProto model;
    sl_.GetSolver().ExportModel(&model);
    RegisterUDFBuilder();
    if (!solver_.LoadModel(model)) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Cannot create the validator solver!";
    }

    /*
     * Now, we want to find all decision variables. We do this via a
     * custom visitor.
     */
    VariableFinder var_finder;
    solver_.Accept(&var_finder);
    VariableFinder::StringVarMap var_map = var_finder.GetVarMap();
    for (StringVector::const_iterator cit = var_names.begin();
            cit != var_names.end(); cit++) {
        VariableFinder::StringVarMap::const_iterator var = var_map.find(*cit);
        if (var == var_map.end()) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Cannot find a variable in the duplicate solver: name="
                    << *cit;
        }
        /*
         *  const_cast is perfectly fine. Vars just come from the visitor,
         *  which does not change them. But these vars are supposed to be
         *  changed during the search anyway, albeit later.
         */
        search_vars_prototype_.Add(const_cast<IntVar *>(var->second));
    }
}

void Validator::AddSolution(const Assignment &sol) {
    // All assignments have the same structure
    Assignment *asgn = new Assignment(&search_vars_prototype_);

    Assignment::IntContainer &to_vars = asgn->MutableIntVarContainer();
    const Assignment::IntContainer &from_vars = sol.IntVarContainer();

    for (size_t i = 0; i < to_vars.Size(); i++) {
        const IntVarElement &from_elem = from_vars.Element(i);
        IntVarElement &to_elem = to_vars.MutableElement(i);
        to_elem.SetRange(from_elem.Min(), from_elem.Max());
        // The only thing we ignore here is "Activated", which is not needed
    }

    to_validate_mtx_.lock();
    to_validate_.push_back(AssignmentPtr(asgn));
    to_validate_mtx_.unlock();
    validate_cond_.notify_one();
}

void Validator::operator()() {
    DecisionBuilder *db = solver_.RevAlloc(new RestoreAssignmentBuilder(*this));
    solver_status_ = solver_.Solve(db, &collector_);
}

IntExpr* Validator::UDFBuilder(std::string name, CPModelLoader* const builder,
        const CPIntegerExpressionProto& proto) {

    // should have IntVar parameters protobuffed
    std::vector<IntVar*> vars;
    if (!builder->ScanArguments(ModelVisitor::kVarsArgument, proto, &vars)) {
        return NULL;
    }

    UDFFunctionCreator udf_creator = sl_.GetRegisteredUDF(name);
    if (!udf_creator) {
        return NULL;
    }

    return udf_creator(&solver_, adapter_, vars);
}

void Validator::RegisterUDFBuilder() {
    StringSet udfs = sl_.GetAllUsedUDFs();
    for (StringSet::const_iterator cit = udfs.begin(); cit != udfs.end();
            cit++) {
        solver_.RegisterBuilder(*cit, NewPermanentCallback(this,
                &Validator::UDFBuilder, std::string(*cit)));
        /*
         * std::string() is used to avoid references inside the callback.
         * Maybe we do not need it, but my knowledge of templates is vague
         * to decide for sure.
         */
    }
}

AssignmentPtrVector *Validator::GetNextAssignments() {
    /*
     * For now, a simple policy: wait until we get an assignmeny to
     * validate and then check it.
     */
    boost::unique_lock<boost::mutex> validate_lock(to_validate_mtx_);
    while (to_validate_.empty() && !terminate_) {
        validate_cond_.wait(validate_lock);
    }

    if (terminate_) {
        return NULL;
    }

    AssignmentPtrVector res = new AssignmentPtrVector;
    res.reserve(to_validate_.size());
    res.insert(res.end(), to_validate_.begin(), to_validate_.end());
    to_validate_.clear();

    return res;
}

} /* namespace searchlight */

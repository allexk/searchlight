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
 * @file ortools_inc.h
 * This is common file for other sources using Google's or-too;s. It contains
 * common definitions and includes from or-tools.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_ORTOOLS_INC_H_
#define SEARCHLIGHT_ORTOOLS_INC_H_

#include <constraint_solver/constraint_solver.h>
#include <constraint_solver/constraint_solveri.h>

#include "base.h"

namespace searchlight {

using operations_research::Solver;

using operations_research::DecisionBuilder;
using operations_research::Decision;
using operations_research::Assignment;
using operations_research::ModelVisitor;
using operations_research::Demon;
using operations_research::SearchMonitor;
using operations_research::SolutionCollector;
using operations_research::SearchLimit;

using operations_research::IntVar;
using operations_research::IntVarIterator;
using operations_research::InitAndGetValues;
using operations_research::IntVarElement;
using operations_research::IntervalVar;
using operations_research::SequenceVar;

using operations_research::BaseIntExpr;
using operations_research::IntExpr;

using operations_research::Constraint;

using operations_research::StringPrintf;
using operations_research::ACMRandom;

/**
 * A vector of Int Variables.
 */
typedef std::vector<IntVar *> IntVarVector;

/**
 * Pointer to a full assignment.
 */
typedef std::unique_ptr<Assignment> AssignmentPtr;


/**
 * Converts a lite assignment to a full Assignment.
 *
 * @param asgn full assignment
 * @param lite_asgn lite assignment
 */
void LiteToFullAssignment(Assignment &asgn,
        const LiteVarAssignment &lite_asgn);

/**
 * Converts a full Assignment to a lite assignment.
 *
 * Elements are added starting from position pos.
 *
 * @param asgn full assignment
 * @param pos starting position for putting elements
 * @param lite_asgn lite assignment
 */
void FullAssignmentToLite(const Assignment &asgn, size_t pos,
        LiteVarAssignment &lite_asgn);

} /* namespace searchlight */

#endif /* SEARCHLIGHT_ORTOOLS_INC_H_ */

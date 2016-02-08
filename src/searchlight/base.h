/*
 * Copyright 2013, Brown University, Providence, RI.
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
 * @file base.h
 * This file contains some basic definitions shared by multiple sources.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_BASE_H_
#define SEARCHLIGHT_BASE_H_

#include <stdint.h>

#include <map>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>

#include <log4cxx/log4cxx.h>
#include <boost/property_tree/ptree.hpp>

namespace searchlight {

/**
 * Property tree used to store configuration options.
 */
using SearchlightConfig = boost::property_tree::ptree;

/**
 * A general integer interval.
 */
typedef std::pair<int64_t, int64_t> IntInterval;

/**
 * Vector of int64 values.
 */
typedef std::vector<int64_t> Int64Vector;

/**
 * A general double interval.
 */
typedef std::pair<double, double> DoubleInterval;

/**
 * A vector of strings.
 */
typedef std::vector<std::string> StringVector;

/**
 * A set of strings.
 */
typedef std::set<std::string> StringSet;

/**
 * Vector of int64
 */
typedef std::vector<int64_t> Int64Vector;

/**
 * Vector of typical sizes.
 */
typedef std::vector<size_t> SizeVector;

/**
 * Vector of double values.
 */
typedef std::vector<double> DoubleVector;

/**
 * An assignment to a variable.
 *
 * It is called "lite" to contrast the full-feature or-tools Assignment,
 * which contains more features when we need most of the time.
 *
 * maxs_ part might be empty, if the assignment is complete, i.e., not a range,
 * but a single value for each variable.
 */
struct LiteVarAssignment {
    /**
     * Min/max values.
     */
    Int64Vector mins_, maxs_;

    /**
     * Clears the entire assignment.
     */
    void clear() {
        mins_.clear();
        maxs_.clear();
    }
};

/**
 *  A vector of lite assignments.
 */
typedef std::vector<LiteVarAssignment> LiteAssignmentVector;

/**
 *  Candidate assignment, including metadata.
 */
struct CandidateAssignment {
    /**
     * The assignment.
     */
    LiteVarAssignment var_asgn_;

    /**
     * Relaxed constraints (const. id, left param., right param.)
     */
    Int64Vector relaxed_constrs_;

    /**
     * Id >= 0, the candidate is remote; -1, needs simulation; -2, local.
     *
     * When we send candidates for forward, forw_id_ >=0
     */
    int64_t forw_id_;
};

/**
 * Vector of candidate assignments.
 */
using CandidateVector = std::vector<CandidateAssignment>;

} /* namespace searchlight */
#endif /* SEARCHLIGHT_BASE_H_ */

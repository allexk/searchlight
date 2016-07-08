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
#include <boost/tokenizer.hpp>

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
 * Bool vector.
 */
typedef std::vector<bool> BoolVector;

/**
 * Check for domination.
 *
 * @param v1 first vector
 * @param v2 second vector
 * @param maxim true, if the component is ascending; false, otherwise
 * @return 1, v1 dominates v2; -1, is dominated; 0 cannot say
 */
int RelateVectors(const Int64Vector &v1, const Int64Vector &v2,
                  const std::vector<bool> &maxim);

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
     * Check if this is a range.
     *
     * @return true, if range; false, otherwise
     */
    bool IsRange() const {
        return !maxs_.empty();
    }

    /**
     * Return size of the assignment.
     *
     * @return assignment size
     */
    size_t Size() const {
        return mins_.size();
    }

    /**
     * Check for domination.
     *
     * @param v vector to check
     * @param maxim true, if the component is ascending; false, otherwise
     * @return 1, this assignment dominates v; -1, is dominated; 0 cannot say
     */
    int Relate(const Int64Vector &v, const std::vector<bool> &maxim) const {
        assert(v.size() == mins_.size() && mins_.size() == maxim.size());
        if (!IsRange()) {
            return RelateVectors(mins_, v, maxim);
        }
        // Range case
        int res = 0;
        for (size_t i = 0; i < v.size(); ++i) {
            int c_res = 0;
            if (v[i] > maxs_[i] || v[i] < mins_[i]) {
                c_res = mins_[i] > v[i] ? 1 : -1;
                // Invert if we're minimizing
                if (!maxim[i]) {
                    c_res = 0 - c_res;
                }
            } else if (mins_[i] != maxs_[i]) { // Really an interval
                // For ranges, if v[i] is within the interval, we cannot say
                return 0;
            }
            // If expected and component are different, cannot say
            if (res && c_res && res != c_res) {
                return 0;
            }
            if (!res) {
                res = c_res;
            }
        }
        return res;
    }

    /**
     * Clears the entire assignment.
     */
    void clear() {
        mins_.clear();
        maxs_.clear();
    }

    /**
     * Check if the assignment is empty.
     *
     * @return true, if empty
     */
    bool empty() const {
        return mins_.empty();
    }

    /**
     * Adds a value pair to the vector.
     *
     * @param val_min minimum value
     * @param val_max maximum value
     */
    void add_val(int64_t val_min, int64_t val_max) {
        mins_.push_back(val_min);
        maxs_.push_back(val_max);
    }
};
std::ostream &operator<<(std::ostream &os, const LiteVarAssignment &lva);

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

    /**
     * Best relaxation degree possible for the assignment.
     */
    double best_rd_;

    /**
     * Best rank possible.
     *
     * Makes sense only if best_rd is 0.0.
     */
    double best_rank_;

    /**
     * Vector of relaxation constraint values.
     *
     */
    LiteVarAssignment rc_vals_;
};
std::ostream &operator<<(std::ostream &os, const CandidateAssignment &ca);

/**
 * Vector of candidate assignments.
 */
using CandidateVector = std::vector<CandidateAssignment>;

/**
 * Approximate value, returned by an estimator.
 *
 * TODO: there might be some benefits in returning the probability of MAY_NULL
 * in the future.
 */
struct IntervalValue {
    /**
     * Describes the state of the value
     */
    enum State {
        NON_NULL,//!< NON_NULL definitely not NULL
        MAY_NULL,//!< MIN_NULL might be NULL (thus, min is invalid)
        NUL      //!< NUL definitely NULL
    };

    /**
     * Minimum possible
     */
    double min_;

    /**
     * Maximum possible
     */
    double max_;

    /**
     * An approximate value.
     */
    double val_;

    /**
     * Value state
     */
    State state_;

    /**
     * Creates a NULL interval value.
     */
    IntervalValue() : min_(0), max_(0), val_(-1), state_(NUL) {}
};

/**
 * Outputs an IntervalValue into an output stream.
 *
 * @param os the stream to output to
 * @param iv the value to output
 * @return the stream used for the output
 */
std::ostream &operator<<(std::ostream &os, const IntervalValue &iv);

/**
 * A vector of interval values.
 */
typedef std::vector<IntervalValue> IntervalValueVector;

/**
 * String tokenizer.
 *
 * @param str string to tokenize
 * @param seps separator characters
 * @return resulting tokens
 */
StringVector TokenizeString(const std::string &str, const char *seps);

} /* namespace searchlight */
#endif /* SEARCHLIGHT_BASE_H_ */

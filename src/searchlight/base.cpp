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
 * @file base.cpp
 * This file contains some basic implementations.
 *
 * @author Alexander Kalinin
 */

#include "base.h"

namespace searchlight {

int RelateVectors(const Int64Vector &v1, const Int64Vector &v2,
                  const std::vector<bool> &maxim) {
    assert(v1.size() == v2.size() == maxim.size());
    int res = 0;
    for (size_t i = 0; i < v1.size(); ++i) {
        int c_res = 0;
        if (v1[i] != v2[i]) {
            c_res = v1[i] > v2[i] ? 1 : -1;
            // Invert if we're minimizing
            if (!maxim[i]) {
                c_res = 0 - c_res;
            }
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

std::ostream &operator<<(std::ostream &os, const IntervalValue &iv) {
    if (iv.state_ == IntervalValue::NUL) {
        os << "(NULL)";
    } else {
        os << "(Min: " << iv.min_ << ", Max: " << iv.max_ <<
                ", Val: " << iv.val_;
        if (iv.state_ == IntervalValue::MAY_NULL) {
            os << ", null?";
        }
        os << ")";
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const LiteVarAssignment &lva) {
    os << "Mins: (";
    std::copy(lva.mins_.begin(), lva.mins_.end(),
              std::ostream_iterator<int64_t>(os, ","));
    os << "), Maxs:(";
    std::copy(lva.maxs_.begin(), lva.maxs_.end(),
              std::ostream_iterator<int64_t>(os, ","));
    os << ")";
    return os;
}

std::ostream &operator<<(std::ostream &os, const CandidateAssignment &ca) {
    os << ca.var_asgn_;
    if (!ca.relaxed_constrs_.empty()) {
        os << ", RC:(";
        std::copy(ca.relaxed_constrs_.begin(), ca.relaxed_constrs_.end(),
                  std::ostream_iterator<int64_t>(os, ","));
        os << ")";
    }
    os << "FW=" << ca.forw_id_ << ", RD=" << ca.best_rd_;
    return os;
}

StringVector TokenizeString(const std::string &str, const char *seps) {
    StringVector res;
    using TokenSeparator = boost::char_separator<char>;
    using Tokenizer = boost::tokenizer<TokenSeparator>;
    TokenSeparator sep{seps}; // parts are separated by '_'
    Tokenizer tokenizer{str, sep};

    for (auto cit = tokenizer.begin(); cit != tokenizer.end(); ++cit) {
        res.push_back(*cit);
    }
    return res;
}
} /* namespace searchlight */

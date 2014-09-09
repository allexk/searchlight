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
 * @file ortools_inc.cpp
 *
 * This file contains some or-tools related stuff. For example, converting
 * between or-tools and Searchlight structures,
 * like Assignment <-> LiteAssignment.
 *
 * @author Alexander Kalinin
 */

#include "ortools_inc.h"

namespace searchlight {

void LiteToFullAssignment(Assignment &asgn,
        const LiteVarAssignment &lite_asgn) {
    assert(lite_asgn.maxs_.empty() ||
            lite_asgn.mins_.size() == lite_asgn.maxs_.size());

    auto &asgn_vars = *asgn.MutableIntVarContainer();
    const bool lite_is_range = !lite_asgn.maxs_.empty();

    assert(asgn_vars.Size() == lite_asgn.mins_.size());

    for (int i = 0; i < asgn_vars.Size(); i++) {
        auto &var = *asgn_vars.MutableElement(i);
        var.SetMin(lite_asgn.mins_[i]);
        if (lite_is_range) {
            var.SetMax(lite_asgn.maxs_[i]);
        }
    }
}

void FullAssignmentToLite(const Assignment &asgn,
        LiteVarAssignment &lite_asgn) {
    assert(lite_asgn.maxs_.empty() ||
            lite_asgn.mins_.size() == lite_asgn.maxs_.size());

    const auto &asgn_vars = asgn.IntVarContainer();
    lite_asgn.mins_.resize(asgn_vars.Size());
    lite_asgn.maxs_.resize(asgn_vars.Size());

    for (int i = 0; i < asgn_vars.Size(); i++) {
        const auto &var = asgn_vars.Element(i);
        lite_asgn.mins_[i] = var.Min();
        lite_asgn.maxs_[i] = var.Max();
    }
}
} /* namespace searchlight */

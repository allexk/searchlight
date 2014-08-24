/* Copyright 2014, Brown University, Providence, RI.
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
 * @file ortools_model.h
 *
 * This file contains some useful methods for manipulating or-tools' CP models.
 *
 * The first one is model duplication. This is needed by the validator and
 * might be needed to replicate solver to multiple threads (think multi-core
 * search). or-tools itself has everything that is needed, but that stuff
 * is hidden inside its implementation. The corresponding cpp file contains
 * copy-pasted and slightly modified classes from the implementation to make
 * model duplication easier.
 *
 * @author Alexander Kalinin
 * @author or-tools Team
 */

#ifndef SEARCHLIGHT_ORTOOLS_MODEL_H_
#define SEARCHLIGHT_ORTOOLS_MODEL_H_

#include "searchlight.h"

namespace searchlight {

/**
 * Clones model from one solver to another.
 *
 * Model consists of variables and constraints. The model is cloned by first
 * exporting it from one solver into a protobuf, and then importing it into
 * another solver. The destination solver should be empty. Otherwise, conflicts
 * are possible and the result is undefined.
 *
 * @param sl searchlight instance
 * @param from solver to get the model from
 * @param to solver to clone the model to
 * @param adapter adapter for UDFs in the destination solver
 * @return true, if clone is successful; false, otherwise
 */
bool CloneModel(const Searchlight &sl, const Solver &from, Solver &to,
        const AdapterPtr &adapter);
}
#endif /* SEARCHLIGHT_ORTOOLS_MODEL_H_ */

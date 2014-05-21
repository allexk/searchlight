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
 * @file common.h
 *
 * This file contains definitions for some common functions and classes useful
 * for writing searchlight tasks
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_TASKS_COMMON_H_
#define SEARCHLIGHT_TASKS_COMMON_H_

#include "searchlight/base.h"

namespace searchlight {

/**
 * Map describing parameters: key (string) -> value (string)
 */
typedef std::map<std::string, std::string> ParamsMap;

/**
 * Parses parameter string and extracts (key, value) pairs. The string
 * is supposed to be in the format: "ke1=value1,key2=value2,...".
 *
 * @param params parameters string
 * @param pmap map to fill with (key, value) pairs
 */
void ParseParameters(const std::string &params, ParamsMap &pmap);


} /* namespace searchlight */
#endif /* SEARCHLIGHT_TASKS_COMMON_H_ */

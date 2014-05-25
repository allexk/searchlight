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

#include <log4cxx/log4cxx.h>

namespace searchlight {

/**
 * A general integer interval.
 */
typedef std::pair<int64_t, int64_t> IntInterval;

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

} /* namespace searchlight */
#endif /* SEARCHLIGHT_BASE_H_ */

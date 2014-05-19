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
 * @file common.cpp
 *
 * Implementation of common taks function and classes.
 *
 * @author Alexander Kalinin
 */

#include "common.h"

#include <boost/tokenizer.hpp>

namespace searchlight {

void ParseParameters(const std::string &params, ParamsMap &pmap) {
    typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
    tokenizer tokens(params, boost::char_separator<char>(";|,"));
    searchlight::StringVector str_tokens;
    for (tokenizer::iterator it = tokens.begin(); it != tokens.end();
            ++it) {
        size_t eq_pos = it->find('=');
        if (eq_pos != std::string::npos && eq_pos != it->size() - 1) {
            std::string key = it->substr(0, eq_pos);
            std::string val = it->substr(eq_pos + 1);
            pmap[key] = val;
        }
    }
}
} /* namespace searchlight */

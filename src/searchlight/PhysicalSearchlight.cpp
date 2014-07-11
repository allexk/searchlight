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
 * @file PhysicalSearchlight.cpp
 * The implementation of the physical part of the searchlight operator.
 *
 * @see LogicalSearchlight.cpp for the description of the operator itself.
 *
 * @author Alexander Kalinin
 */

#include "searchlight_task.h"

#include <query/Operator.h>
#include <boost/tokenizer.hpp>

namespace scidb
{

/**
 * Implements the physical part of the operator. It instantiates a new
 * SearchlightTask, which calls an appropriate library function. It also
 * returns the resulting array, which is lazily filled by calls from the
 * upstream operators.
 */
class PhysicalSearchlight : public PhysicalOperator
{
public:
    /**
     * Constructs a new physical part of the Searchlight operator.
     *
     * @param logicalName the name of the logical part
     * @param physicalName the name of the physical part
     * @param parameters params (parsed in the logical part)
     * @param schema the schema of the result
     */
    PhysicalSearchlight(const std::string &logicalName,
            const std::string &physicalName, const Parameters &parameters,
            const ArrayDesc &schema) :
        PhysicalOperator(logicalName, physicalName, parameters, schema) {}

    /**
     * Execute the operator and return the output array. Note that we do not
     * run the actual search until the first chunk of the array is requested.
     *
     * @param inputArrays the input array arguments. Data and sample.
     * @param query the query context
     * @return the output array object
     */
    boost::shared_ptr<Array> execute(
            std::vector<boost::shared_ptr<Array>> &inputArrays,
            boost::shared_ptr<Query> query) {

        // tokenize params
        std::string lib_func_param =
                ((boost::shared_ptr<OperatorParamPhysicalExpression> &)
                        _parameters[0])->getExpression()->
                        evaluate().getString();
        typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
        tokenizer tokens(lib_func_param, boost::char_separator<char>(":", 0));
        searchlight::StringVector str_tokens;
        for (tokenizer::iterator it = tokens.begin(); it != tokens.end();
                ++it) {
            str_tokens.push_back(*it);
        }

        // FIXME: should probably check this in the logical part
        if (str_tokens.size() != 3) {
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                    SCIDB_LE_ILLEGAL_OPERATION) <<
                            "Illegal searchlight parameters!";
        }

        // Gather all samples
        ArrayPtrVector samples;
        for (size_t i = 1; i < inputArrays.size(); i++) {
            samples.push_back(inputArrays[i]);
        }
        assert(!samples.empty());

        searchlight::SearchlightTaskPtr sl_task =
                boost::make_shared<searchlight::SearchlightTask>(str_tokens[0],
                        str_tokens[1], str_tokens[2]);
        sl_task->GetSearchlight().RegisterArray(inputArrays[0], samples);
        return ArrayPtr(new searchlight::SearchlightResultsArray(sl_task,
                _schema, query));
    }
};

/*
 * Register the operator.
 */
REGISTER_PHYSICAL_OPERATOR_FACTORY(PhysicalSearchlight, "searchlight",
        "PhysicalSearchlight");

} /* namespace scidb */

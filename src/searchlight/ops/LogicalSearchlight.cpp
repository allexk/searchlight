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
 * @file LogicalSearchlight.cpp
 * The implementation of the logical part of the searchlight operator.
 *
 * @brief The operator: searchlight().
 *
 * @par Synopsis:
 *   <br/>
 *   searchlight(input_array, grid_array [, grid_array]*,
 *      'library:function:task_json_spec').
 *
 * @par Summary:
 *   <br/>
 *   This operator runs a task via Searchlight. It finds the specified function
 *   in the specified library, which is a DLL. The function itself takes
 *   the Searechlight instance and a string of parameters, which it supposed
 *   to parse by itself.
 *   <br/>
 *   The operator outputs an array of results with indexes [0, *]. The right
 *   boundary is undefined since the number of results is unknown beforehand.
 *   The array is output in the stream mode, which allows it to be processed
 *   online.
 *
 * @par Input:
 *   <br> input_array -- any SciDB array
 *   <br> grid_array -- grid arrays for the specified array. The array
 *      must correspond to the following format: dimensions[chunk: 0, *;
 *      sampled_attribute: original attribute index], attributes:[min:double,
 *      max: double, count: uint64, sum: double]. The user can specify several
 *      of them, in which case the first one will be used as the sample, and
 *      others will be used to optimize aggregate windows computation. It is
 *      assumed that array are listed from lower to higher "resolutions".
 *      The name of each must contain suffix "_NxNx...", where N specifies the
 *      step of the grid for each dimension.
 *   <br> 'library:function:task_json_spec'  -- the name of the DLL where the
 *      tasks come from. The name of the task (function), which should take
 *      a ref to Searchlight and a path to the JSON taks specification. It will
 *      be parsed by the task itself.
 *
 * @par Output array:
 *   <br> <
 *   <br>   assignment :string
 *   <br> >
 *   <br> [
 *   <br>   num = 0:*,1 0]
 *   <br> ]
 *
 * @author Alexander Kalinin
 */

#include "../searchlight_task.h"
#include <query/Operator.h>

namespace scidb {

/**
 * The Logical Operator object for searchlight.
 */
class LogicalSearchlight : public LogicalOperator {
public:
    /**
     * Searchlight operator constructor.
     *
     * @param logicalName used internally by scidb
     * @param alias used internally by scidb
     */
    LogicalSearchlight(const std::string& logicalName,
            const std::string& alias) :
        LogicalOperator(logicalName, alias) {

        _properties.tile = true; // enable propagation of tile mode

        ADD_PARAM_INPUT(); // array
        ADD_PARAM_INPUT(); // primary sample
        ADD_PARAM_VARIES();
    }

    /**
     * Returns options for the next parameter in the variable parameter list
     * @param schemas input schemas parsed before
     * @return a vector of possible options
     */
    virtual std::vector<boost::shared_ptr<OperatorParamPlaceholder>>
        nextVaryParamPlaceholder(
                const std::vector<ArrayDesc> &schemas) override {

        std::vector<boost::shared_ptr<OperatorParamPlaceholder>> res;
        if (_parameters.empty()) {
            // aux sample or the task spec
            res.push_back(PARAM_INPUT());
            res.push_back(PARAM_CONSTANT(TID_STRING));
        } else {
            // we already have the task spec (inputs are not added to params)
            res.push_back(END_OF_VARIES_PARAMS());
        }

        return res;
    }

    /**
     * Returns the descriptor of the result.
     *
     * @param schemas all of the schemas of the input arrays
     *  (if the operator accepts any)
     * @param query the query context
     * @return the schema of the output
     */
    virtual ArrayDesc inferSchema(std::vector<ArrayDesc> schemas,
            boost::shared_ptr<Query> query) override {

        return searchlight::SearchlightResultsArray::GetArrayDesc();
    }
};

// Register the operator
REGISTER_LOGICAL_OPERATOR_FACTORY(LogicalSearchlight, "searchlight");

} /* namespace scidb */

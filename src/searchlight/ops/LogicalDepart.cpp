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
 * @file LogicalDepart.cpp
 * The implementation of the logical part of the departitioner operator.
 *
 * @brief The operator: depart().
 *
 * @par Synopsis:
 *   <br/>
 *   depart(input_array, [bool prereplicate]).
 *
 * @par Summary:
 *   <br/>
 *   This operator "de-partitions" the underlying array, effectively
 *   making lazy, on-demand replication. The schema remains the same, and
 *   all chunks are returned, even if those not present at the instance.
 *   The operator determines where the absent chunk should be located and
 *   requests a transfer. The operator is supposed to be used with the
 *   searchlight() operator, and it uses searchlight mesenger for requests and
 *   responses. If the array is already replicated, it does nothing.
 *
 * @par Input:
 *   <br/> input_array -- any SciDB array
 *   <br/> pre_replicate -- if true, then the operator will request replication
 *         from SciDb. Otherwise, it will request HashPartitioned distribution,
 *         which is the SciDb's default for most arrays. The parameter itself
 *         is false by default.
 *
 * @par Output array:
 *   <br/> the same as input array
 *
 * @author Alexander Kalinin
 */

#include "../scidb_inc.h"
#include <query/Operator.h>

namespace scidb {

/**
 * The Logical Operator object for depart.
 */
class LogicalDepart : public LogicalOperator {
public:
    /**
     * Depart operator constructor.
     *
     * @param logicalName used internally by scidb
     * @param alias used internally by scidb
     */
    LogicalDepart(const std::string& logicalName,
            const std::string& alias) :
        LogicalOperator(logicalName, alias) {

        _properties.tile = true; // enable propagation of tile mode
        ADD_PARAM_INPUT(); // array
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
            res.push_back(PARAM_CONSTANT(TID_BOOL));
        }
        res.push_back(END_OF_VARIES_PARAMS());

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
        assert(schemas.size() == 1);
        if (schemas[0].getSize() == scidb::INFINITE_LENGTH) {
            throw USER_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Array for departitioning must be finite!";
        }

        return schemas[0];
    }
};

// Register the operator
REGISTER_LOGICAL_OPERATOR_FACTORY(LogicalDepart, "depart");

} /* namespace scidb */

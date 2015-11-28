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
 *   depart(input_array, [lb+, rb+], [bool prereplicate]).
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
 *   Additionally the user can specify the array window restriction to convey
 *   additional information about the part of the array accessed to
 *   Searchlight. This information will be conveyed via the array's schema.
 *
 * @par Input:
 *   <br/> input_array -- any SciDB array
 *   <br/> lb+, rb+ -- window coordinates
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
        assert(schemas.size() == 1);
        const size_t params_num = _parameters.size();
        if (params_num == 0) {
            res.push_back(PARAM_CONSTANT(TID_BOOL));
            res.push_back(PARAM_CONSTANT(TID_INT64));
            res.push_back(END_OF_VARIES_PARAMS());
        } else if (params_num == 1) {
            assert(_parameters[0]->getParamType() == PARAM_LOGICAL_EXPRESSION);
            const auto &param =
                    boost::static_pointer_cast<OperatorParamLogicalExpression>
                        (_parameters[0]);
            if (param->getExpectedType().typeId() == TID_BOOL) {
                res.push_back(END_OF_VARIES_PARAMS());
            } else {
                assert(param->getExpectedType().typeId() == TID_INT64);
                res.push_back(PARAM_CONSTANT(TID_INT64));
            }
        } else {
            const size_t dims_num = schemas[0].getDimensions().size();
            if (params_num < 2 * dims_num) {
                res.push_back(PARAM_CONSTANT(TID_INT64));
            } else if (params_num == 2 * dims_num) {
                res.push_back(PARAM_CONSTANT(TID_BOOL));
                res.push_back(END_OF_VARIES_PARAMS());
            } else {
                res.push_back(END_OF_VARIES_PARAMS());
            }
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
        assert(schemas.size() == 1);
        ArrayDesc &array_desc = schemas[0];

       /*
        * Trim the descriptor. We need this for properly functioning
        * Searchlight code. We probably care only about the finite
        * lower boundary.
        *
        * FIXME: Check if we really need this, especially in the sampler.
        */
       array_desc.trim();

        if (_parameters.size() > 1) {
            // Copy should work, but the descriptor pointer will be incorrect.
            // It will change when we create a new descriptor, though.
            const Dimensions &dims = array_desc.getDimensions();
            Dimensions new_dims;
            for (size_t i = 0; i < dims.size(); i++) {
                const auto &param_low =
                    boost::static_pointer_cast<OperatorParamLogicalExpression>
                            (_parameters[i]);
                const auto &param_high =
                    boost::static_pointer_cast<OperatorParamLogicalExpression>
                            (_parameters[i + dims.size()]);
                Coordinate low = evaluate(param_low->getExpression(),
                        query, TID_INT64).getInt64();
                Coordinate high = evaluate(param_high->getExpression(),
                        query, TID_INT64).getInt64();

                // First, we want to check the boundaries. Note, we cannot
                // have a infinite length array here.
                const DimensionDesc &dim = dims[i];
                if (low < dim.getStartMin()) {
                    low = dim.getStartMin();
                }
                if (high > dim.getEndMax()) {
                    high = dim.getEndMax();
                }

                // Then, we want to align the coordinates with the chunks
                Coordinate offs = low - dim.getStartMin();
                if (offs % dim.getChunkInterval()) {
                    low -= offs % dim.getChunkInterval();
                }
                offs = high - dim.getStartMin();
                if ((offs + 1) % dim.getChunkInterval()) {
                   high -= offs % dim.getChunkInterval();
                   high += (dim.getChunkInterval() - 1);
                }

                // Correctness
                if (low > high) {
                    throw USER_EXCEPTION(
                            SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                                    << "Incorrect coordinates!";
                }
                new_dims.emplace_back(dim.getBaseName(),
                        dim.getNamesAndAliases(), low, low, high, high,
                        dim.getChunkInterval(), dim.getChunkOverlap());
            }

            return ArrayDesc(array_desc.getName(), array_desc.getAttributes(),
                    new_dims, array_desc.getFlags());
        } else {
            return array_desc;
        }
    }

private:
    // Window, if specified
    Coordinates low_, high_;
};

// Register the operator
REGISTER_LOGICAL_OPERATOR_FACTORY(LogicalDepart, "depart");

} /* namespace scidb */

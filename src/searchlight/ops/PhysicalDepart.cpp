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
 * @file PhysicalDepart.cpp
 * The implementation of the physical part of the depart operator.
 *
 * @see LogicalDepart.cpp for the description of the operator itself.
 *
 * @author Alexander Kalinin
 */

#include "../scidb_inc.h"
#include <query/Operator.h>

#include "DepartArray.h"

namespace scidb
{

/**
 * Implements the physical part of the operator. It creates a new array that
 * either returns local chunks or requests, materializes and returns remote
 * chunks.
 */
class PhysicalDepart : public PhysicalOperator
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
    PhysicalDepart(const std::string &logicalName,
            const std::string &physicalName, const Parameters &parameters,
            const ArrayDesc &schema) :
        PhysicalOperator(logicalName, physicalName, parameters, schema) {

        bool require_replicate = false;
        if (_parameters.size() > 0) {
                require_replicate =
                    ((boost::shared_ptr<OperatorParamPhysicalExpression> &)
                        _parameters[0])->getExpression()->
                        evaluate().getBool();
        }
        if (require_replicate) {
            input_distr_ = ArrayDistribution(psReplication);
        } else {
            input_distr_ = ArrayDistribution(psHashPartitioned);
        }
    }

    /**
      * [Optimizer API] Determine if operator changes result chunk distribution.
      * @param  input_schemas schemas for all inputs of the operator
      * @return true, if change the distribution; false, otherwise
      */
    virtual bool changesDistribution(
            const std::vector<ArrayDesc> &input_schemas) const override {
        /*
         * We want to change the distribution, since Habilis tend to insist
         * on psHashPartitioned at times, especially when an input is
         * replicated. While this makes sense sometimes, and, actually, we are
         * going to output psReplicated, we want to avoid Habilis to be
         * messing with our output. For example, by inserting reduce_distro().
         *
         * That is why we change our distribution to psUndefined. This should
         * be fine with the searchlight operator and many others. However,
         * if the parent
         * operator requests specific distribution, this will result in
         * inserting parent sg() by Habilis and unnecessary redistribution.
         * This is hard (impossible?) to avoid right now.
         */
        return true;
    }

    /**
     *  [Optimizer API] Determine the distribution of operator output.
     *  @param source_distrs input distributions
     *  @param source_schemas input schemas
     *  @return distribution of the output
     */
    virtual ArrayDistribution getOutputDistribution(
            const std::vector<ArrayDistribution> &source_distrs,
            const std::vector<ArrayDesc> &source_schemas) const override {
        // see the comment for changeDistribution()
        return ArrayDistribution(psUndefined);
    }

    /**
     * Specifies the distribution requirement from the operator to its
     * children.
     *
     * @param source_schemas input schemas
     * @return the required distribution
     */
    virtual DistributionRequirement getDistributionRequirement(
            const std::vector<ArrayDesc> &source_schemas) const override {
        /*
         * What we want to do here is to check if the user wants replication.
         * If she does, enforce it via SciDb's Habilis. Otherwise, ask for
         * HashPartitioned, which is the best alternative, probably and suits
         * most of the arrays flying around here.
         */
        std::vector<ArrayDistribution> reqs{input_distr_};
        return DistributionRequirement(
                DistributionRequirement::SpecificAnyOrder,
                reqs);
    }

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
        if (input_distr_.getPartitioningSchema() == psReplication) {
            return inputArrays[0];
        } else {
            boost::shared_ptr<Query> query = Query::getValidQueryPtr(_query);
            searchlight::SearchlightMessenger::getInstance()->
                    RegisterQuery(query);
            return boost::shared_ptr<Array>(
                    new searchlight::DepartArray(_schema, inputArrays[0],
                            input_distr_, query));
        }
    }

private:
    // The distribution of the input
    ArrayDistribution input_distr_;
};

/*
 * Register the operator.
 */
REGISTER_PHYSICAL_OPERATOR_FACTORY(PhysicalDepart, "depart",
        "PhysicalDepart");

} /* namespace scidb */

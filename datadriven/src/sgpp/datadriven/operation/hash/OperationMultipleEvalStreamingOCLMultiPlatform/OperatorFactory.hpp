// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Factory method for creating a new instance of this variant of OperationMultipleEval.
 * This factory method configures the resulting object by using the parameters of the provided
 * configuration.
 * If no parameters are provided, the default parameter values are used.
 * If a configuration is provided, but some entries were not set, this class adds the missing
 * entries with their default values.
 * This class is the non-templated entry point for the templated inner objects.
 * Templates are used to implement different floating point precision.
 *
 * @see OperationMultiEvalStreamingOCLMultiPlatform
 * @param grid The sparse grid to evaluate
 * @param dataset The datapoints to evaluate
 * @param configuration Configuration that may contain a parameter object for configuration details
 */
base::OperationMultipleEval* createStreamingOCLMultiPlatformConfigured(
    base::Grid& grid, base::DataMatrix& dataset,
    sgpp::datadriven::OperationMultipleEvalConfiguration& configuration);

}  // namespace datadriven
}  // namespace sgpp

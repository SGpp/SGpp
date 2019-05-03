// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef QUADRATUREOPFACTORY_HPP_
#define QUADRATUREOPFACTORY_HPP_

/*
 * This file contains factory methods for operations.
 */

#include <sgpp/quadrature/operation/hash/OperationQuadratureMCAdvanced.hpp>
#include <sgpp/globaldef.hpp>

#include <random>

namespace sgpp {
namespace op_factory {

/**
 * Creates an OperationQuadratureMCAdvanced.
 *
 * @param grid Reference to the grid object
 * @param numberOfSamples Number of Monte Carlo samples
 * @param seed Custom seed (defaults to default seed of mt19937_64)
 */
quadrature::OperationQuadratureMCAdvanced* createOperationQuadratureMCAdvanced(
    base::Grid& grid, size_t numberOfSamples, std::uint64_t seed = std::mt19937_64::default_seed);

}  // namespace op_factory
}  // namespace sgpp

#endif /* QUADRATUREOPFACTORY_HPP_ */

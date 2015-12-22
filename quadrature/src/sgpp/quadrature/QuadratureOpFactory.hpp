// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef QUADRATUREOPFACTORY_HPP_
#define QUADRATUREOPFACTORY_HPP_

/*
 * This file contains factory methods for operations.
 */

#include <random>
#include <sgpp/quadrature/operation/hash/OperationQuadratureMCAdvanced.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace op_factory {

    /**
     * Creates an OperationQuadratureMCAdvanced, specifying a grid
     * object and the number of samples to use.
     *
     * @param grid
     * @param numberOfSamples
     */
    quadrature::OperationQuadratureMCAdvanced* createOperationQuadratureMCAdvanced(base::Grid& grid, size_t numberOfSamples,
        std::uint64_t seed = std::mt19937_64::default_seed);

  } /* namespace op_factory */
} /* namespace SGPP */

#endif /* QUADRATUREOPFACTORY_HPP_ */

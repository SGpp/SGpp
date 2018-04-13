// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/quadrature/QuadratureOpFactory.hpp>
#include <sgpp/globaldef.hpp>
#include <random>

namespace sgpp {
namespace op_factory {

quadrature::OperationQuadratureMCAdvanced* createOperationQuadratureMCAdvanced(
    base::Grid& grid, size_t numberOfSamples, std::uint64_t seed) {
  return new quadrature::OperationQuadratureMCAdvanced(grid, numberOfSamples, seed);
}

}  // namespace op_factory
}  // namespace sgpp

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OperationMakePositiveInterpolationAlgorithm.hpp"
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <vector>

namespace sgpp {
namespace base {

OperationMakePositiveInterpolationAlgorithm::OperationMakePositiveInterpolationAlgorithm() {}

OperationMakePositiveInterpolationAlgorithm::~OperationMakePositiveInterpolationAlgorithm() {}

// -------------------------------------------------------------------------------------------

void OperationMakePositiveSetToZero::computeHierarchicalCoefficients(
    base::Grid& grid, base::DataVector& alpha, std::vector<size_t> addedGridPoints) {
  base::HashGridStorage& gridStorage = grid.getStorage();

  // compute the nodal values of the newly added grid points and subtract the
  // nodal value from the hierarchical coefficient. This sets the function value at that point
  // to be zero
  auto opEval = op_factory::createOperationEval(grid);
  base::DataVector x(gridStorage.getDimension());

  for (auto& i : addedGridPoints) {
    gridStorage.getCoordinates(gridStorage.getPoint(i), x);
    alpha[i] -= opEval->eval(alpha, x);
  }
}

} /* namespace base */
} /* namespace sgpp */

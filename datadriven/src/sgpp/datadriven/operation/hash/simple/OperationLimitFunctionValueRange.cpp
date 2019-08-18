// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationLimitFunctionValueRange.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationMakePositive.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <vector>
#include <limits>

namespace sgpp {
namespace datadriven {

OperationLimitFunctionValueRange::OperationLimitFunctionValueRange(
    datadriven::MakePositiveCandidateSearchAlgorithm candidateSearch,
    datadriven::MakePositiveInterpolationAlgorithm interpolationAlgorithm, bool verbose,
    sgpp::base::ScalarFunction* f)
    : verbose(verbose) {
  opPositive.reset(op_factory::createOperationMakePositive(candidateSearch, interpolationAlgorithm,
                                                           true, verbose, f));
}

OperationLimitFunctionValueRange::~OperationLimitFunctionValueRange() {}

void OperationLimitFunctionValueRange::prepareForLowerLimitation(base::Grid& grid,
                                                                 base::DataVector& alpha,
                                                                 double ylower) {
  // f(x) >= ylower => g(x) := f(x) - ylower >= 0
  addConst(grid, alpha, 1.0, -ylower);
}

void OperationLimitFunctionValueRange::inverseFromLowerLimitation(base::Grid& grid,
                                                                  base::DataVector& alpha,
                                                                  double ylower) {
  // f(x) = g(x) + ylower
  addConst(grid, alpha, 1.0, ylower);
}

void OperationLimitFunctionValueRange::prepareForUpperLimitation(base::Grid& grid,
                                                                 base::DataVector& alpha,
                                                                 double yupper) {
  // f(x) <= yupper => -f(x) >= -yupper => g(x) := -f(x) + yupper >= 0
  addConst(grid, alpha, -1.0, yupper);
}

void OperationLimitFunctionValueRange::inverseFromUpperLimitation(base::Grid& grid,
                                                                  base::DataVector& alpha,
                                                                  double yupper) {
  // f(x) = -g(x) + yupper
  addConst(grid, alpha, -1.0, yupper);
}

void OperationLimitFunctionValueRange::doLowerLimitation(base::Grid& grid, base::DataVector& alpha,
                                                         double ylower, bool limitNodalValues) {
  // limit the function values from below
  prepareForLowerLimitation(grid, alpha, ylower);
  opPositive->makePositive(grid, alpha, limitNodalValues);
  inverseFromLowerLimitation(grid, alpha, ylower);
}

void OperationLimitFunctionValueRange::doUpperLimitation(base::Grid& grid, base::DataVector& alpha,
                                                         double yupper, bool limitNodalValues) {
  // limit the function from above
  prepareForUpperLimitation(grid, alpha, yupper);
  opPositive->makePositive(grid, alpha, limitNodalValues);
  inverseFromUpperLimitation(grid, alpha, yupper);
}

void OperationLimitFunctionValueRange::doLimitation(base::Grid& grid, base::DataVector& alpha,
                                                    double ylower, double yupper,
                                                    bool limitNodalValues) {
  // limit the nodal values at the grid points
  if (limitNodalValues) {
    // limit the nodal values to the lower range
    prepareForLowerLimitation(grid, alpha, ylower);
    opPositive->makeCurrentNodalValuesPositive(grid, alpha);
    inverseFromLowerLimitation(grid, alpha, ylower);

    // limit the nodal values to the upper range
    prepareForUpperLimitation(grid, alpha, yupper);
    opPositive->makeCurrentNodalValuesPositive(grid, alpha);
    inverseFromUpperLimitation(grid, alpha, yupper);
  }

  size_t oldGridSize = 0;
  while (oldGridSize < grid.getSize()) {
    // store the current grid size
    oldGridSize = grid.getSize();
    // limit the function values from below
    doLowerLimitation(grid, alpha, ylower, false);
    // limit the function values from above
    doUpperLimitation(grid, alpha, yupper, false);
  }
}

void OperationLimitFunctionValueRange::addConst(base::Grid& grid, base::DataVector& alpha, double c,
                                                double y) {
  // dehierarchize the values
  auto opHier = op_factory::createOperationHierarchisation(grid);
  opHier->doDehierarchisation(alpha);
  // add the constant
  for (size_t i = 0; i < alpha.getSize(); i++) {
    alpha[i] = c * alpha[i] + y;
  }
  // hierarchize the result
  opHier->doHierarchisation(alpha);
}

std::vector<size_t>& OperationLimitFunctionValueRange::getAddedGridPoints() {
  return opPositive->getAddedGridPoints();
}

std::vector<size_t>& OperationLimitFunctionValueRange::getAddedGridPointsForRangeLimitation() {
  return opPositive->getAddedGridPointsForPositivity();
}

size_t OperationLimitFunctionValueRange::numAddedGridPoints() {
  return opPositive->numAddedGridPoints();
}

size_t OperationLimitFunctionValueRange::numAddedGridPointsForRangeLimitation() {
  return opPositive->numAddedGridPointsForPositivity();
}

} /* namespace datadriven */
} /* namespace sgpp */

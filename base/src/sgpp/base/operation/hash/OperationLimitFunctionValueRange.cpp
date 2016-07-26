// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationLimitFunctionValueRange.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <limits>

namespace sgpp {
namespace base {

OperationLimitFunctionValueRange::OperationLimitFunctionValueRange(
    base::Grid& grid, base::MakePositiveCandidateSearchAlgorithm candidateSearch,
    base::MakePositiveInterpolationAlgorithm interpolationAlgorithm, bool verbose)
    : grid(grid),
      candidateSearch(candidateSearch),
      interpolationAlgorithm(interpolationAlgorithm),
      verbose(verbose) {}

OperationLimitFunctionValueRange::~OperationLimitFunctionValueRange() {}

void OperationLimitFunctionValueRange::doLowerLimitation(base::Grid*& newGrid,
                                                         base::DataVector& newAlpha,
                                                         double ylower) {
  // limit the function values from below
  // f(x) >= ylower => g(x) := f(x) - ylower >= 0
  addConst(newAlpha, 1.0, -ylower);
  op_factory::createOperationMakePositive(grid, candidateSearch, interpolationAlgorithm, verbose)
      ->makePositive(newGrid, newAlpha);
  // f(x) = g(x) + ylower
  addConst(newAlpha, 1.0, ylower);
}

void OperationLimitFunctionValueRange::doUpperLimitation(base::Grid*& newGrid,
                                                         base::DataVector& newAlpha,
                                                         double yupper) {
  // limit the function from above
  // f(x) <= yupper => -f(x) >= -yupper => g(x) := -f(x) + yupper >= 0
  addConst(newAlpha, -1.0, yupper);
  op_factory::createOperationMakePositive(grid, candidateSearch, interpolationAlgorithm, verbose)
      ->makePositive(newGrid, newAlpha);
  // f(x) = -g(x) + yupper
  addConst(newAlpha, -1.0, yupper);
}

void OperationLimitFunctionValueRange::doLimitation(base::Grid*& newGrid,
                                                    base::DataVector& newAlpha, double ylower,
                                                    double yupper) {
  doLowerLimitation(newGrid, newAlpha, ylower);
  doUpperLimitation(newGrid, newAlpha, yupper);
}

void OperationLimitFunctionValueRange::addConst(base::DataVector& alpha, double c, double y) {
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
} /* namespace base */
} /* namespace sgpp */

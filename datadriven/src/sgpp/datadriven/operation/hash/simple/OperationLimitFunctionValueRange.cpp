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
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <limits>

namespace sgpp {
namespace datadriven {

OperationLimitFunctionValueRange::OperationLimitFunctionValueRange(
    base::Grid& grid, datadriven::MakePositiveCandidateSearchAlgorithm candidateSearch,
    datadriven::MakePositiveInterpolationAlgorithm interpolationAlgorithm,
    bool generateConsistentGrid, bool verbose)
    : grid(grid),
      candidateSearch(candidateSearch),
      interpolationAlgorithm(interpolationAlgorithm),
      generateConsistentGrid(generateConsistentGrid),
      verbose(verbose) {}

OperationLimitFunctionValueRange::~OperationLimitFunctionValueRange() {}

void OperationLimitFunctionValueRange::doLowerLimitation(base::Grid*& newGrid,
                                                         base::DataVector& newAlpha, double ylower,
                                                         bool resetGrid) {
  // limit the function values from below
  // f(x) >= ylower => g(x) := f(x) - ylower >= 0
  if (resetGrid) {
    newGrid = grid.clone();
  }

  addConst(*newGrid, newAlpha, 1.0, -ylower);
  std::unique_ptr<datadriven::OperationMakePositive> opPositive(
      op_factory::createOperationMakePositive(*newGrid, candidateSearch, interpolationAlgorithm,
                                              generateConsistentGrid, verbose));
  opPositive->makePositive(newGrid, newAlpha, false);
  // f(x) = g(x) + ylower
  addConst(*newGrid, newAlpha, 1.0, ylower);
}

void OperationLimitFunctionValueRange::doUpperLimitation(base::Grid*& newGrid,
                                                         base::DataVector& newAlpha, double yupper,
                                                         bool resetGrid) {
  // limit the function from above
  // f(x) <= yupper => -f(x) >= -yupper => g(x) := -f(x) + yupper >= 0
  if (resetGrid) {
    newGrid = grid.clone();
  }

  addConst(*newGrid, newAlpha, -1.0, yupper);
  std::unique_ptr<datadriven::OperationMakePositive> opPositive(
      op_factory::createOperationMakePositive(*newGrid, candidateSearch, interpolationAlgorithm,
                                              generateConsistentGrid, verbose));
  opPositive->makePositive(newGrid, newAlpha, false);
  // f(x) = -g(x) + yupper
  addConst(*newGrid, newAlpha, -1.0, yupper);
}

void OperationLimitFunctionValueRange::doLimitation(base::Grid*& newGrid,
                                                    base::DataVector& newAlpha, double ylower,
                                                    double yupper) {
  size_t newGridSize = grid.getSize();
  size_t oldGridSize = 0;
  size_t iteration = 0;
  while (oldGridSize < newGridSize) {
    oldGridSize = newGridSize;
    doLowerLimitation(newGrid, newAlpha, ylower, iteration == 0);
    doUpperLimitation(newGrid, newAlpha, yupper, false);
    newGridSize = newGrid->getSize();
    iteration++;
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

} /* namespace datadriven */
} /* namespace sgpp */

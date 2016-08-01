// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationLimitFunctionValueRange.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/exception/application_exception.hpp>

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
                                                         base::DataVector& newAlpha, double ylower,
                                                         bool resetGrid) {
  // limit the function values from below
  // f(x) >= ylower => g(x) := f(x) - ylower >= 0
  if (resetGrid) {
    copyGrid(grid, newGrid);
    resetGrid = false;
  }

  addConst(*newGrid, newAlpha, 1.0, -ylower);
  op_factory::createOperationMakePositive(*newGrid, candidateSearch, interpolationAlgorithm,
                                          verbose)
      ->makePositive(newGrid, newAlpha, false);
  // f(x) = g(x) + ylower
  addConst(*newGrid, newAlpha, 1.0, ylower);
}

void OperationLimitFunctionValueRange::doUpperLimitation(base::Grid*& newGrid,
                                                         base::DataVector& newAlpha, double yupper,
                                                         bool resetGrid) {
  // limit the function from above
  // f(x) <= yupper => -f(x) >= -yupper => g(x) := -f(x) + yupper >= 0
  if (resetGrid) {
    copyGrid(grid, newGrid);
  }

  addConst(*newGrid, newAlpha, -1.0, yupper);
  op_factory::createOperationMakePositive(*newGrid, candidateSearch, interpolationAlgorithm,
                                          verbose)
      ->makePositive(newGrid, newAlpha, false);
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

void OperationLimitFunctionValueRange::copyGrid(base::Grid& grid, base::Grid*& newGrid) {
  // create grid of dimensions d - 1 of the same type
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto numDims = gridStorage.getDimension();

  newGrid = base::Grid::createLinearGrid(numDims).release();
  base::HashGridStorage& newGridStorage = newGrid->getStorage();

  // run through grid g and add points to mg
  base::GridPoint gridPoint(numDims);

  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    newGridStorage.insert(gridStorage.getPoint(i));
  }
  newGridStorage.recalcLeafProperty();
}

} /* namespace base */
} /* namespace sgpp */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationMakePositive.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationMakePositiveCandidateSetAlgorithm.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationMakePositiveInterpolationAlgorithm.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>

#include <vector>
#include "../../../../../../../base/src/sgpp/base/function/scalar/ScalarFunction.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace datadriven {

OperationMakePositive::OperationMakePositive(
    MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm,
    MakePositiveInterpolationAlgorithm interpolationAlgorithm, bool generateConsistentGrid,
    bool verbose, sgpp::base::ScalarFunction* f)
    : maxLevel(0),
      minimumLevelSum(0),
      maximumLevelSum(0),
      generateConsistentGrid(generateConsistentGrid),
      candidateSearchAlgorithm(candidateSearchAlgorithm),
      interpolationAlgorithm(interpolationAlgorithm),
      f(f),
      verbose(verbose) {}

OperationMakePositive::~OperationMakePositive() {}

void OperationMakePositive::initialize(base::Grid& grid, base::DataVector& alpha) {
  // set range of level sums for candidate search
  base::HashGridStorage& gridStorage = grid.getStorage();
  size_t numDims = gridStorage.getDimension();

  // init max level for candidate search algorithm
  size_t candidateSearchMaxLevel = gridStorage.getMaxLevel();
  maxLevel = candidateSearchMaxLevel;
  minimumLevelSum = numDims;
  maximumLevelSum = maxLevel * numDims;

  // set the candidate search algorithm
  switch (candidateSearchAlgorithm) {
    case MakePositiveCandidateSearchAlgorithm::FullGrid:
      candidateSearch = std::make_shared<datadriven::OperationMakePositiveLoadFullGridCandidates>(
          candidateSearchMaxLevel);
      break;
    case MakePositiveCandidateSearchAlgorithm::Intersections:
      candidateSearch =
          std::make_shared<datadriven::OperationMakePositiveFindIntersectionCandidates>(
              candidateSearchMaxLevel);
      break;
    case MakePositiveCandidateSearchAlgorithm::HybridFullIntersections:
      candidateSearch =
          std::make_shared<datadriven::OperationMakePositiveHybridFindIntersectionCandidates>(
              candidateSearchMaxLevel);
      break;
    case MakePositiveCandidateSearchAlgorithm::IntersectionsJoin:
      candidateSearch =
          std::make_shared<datadriven::OperationMakePositiveFindIntersectionCandidatesJoin>(
              candidateSearchMaxLevel);
      break;
    default:
      candidateSearch = std::make_shared<datadriven::OperationMakePositiveLoadFullGridCandidates>(
          candidateSearchMaxLevel);
      break;
  }
  candidateSearch->setVerbose(verbose);

  // set the interpolation algorithm
  switch (interpolationAlgorithm) {
    case MakePositiveInterpolationAlgorithm::SetToZero:
      interpolation = std::make_shared<datadriven::OperationMakePositiveSetToZero>();
      break;
    case MakePositiveInterpolationAlgorithm::InterpolateExp:
      interpolation = std::make_shared<datadriven::OperationMakePositiveInterpolateExp>();
      break;
    case MakePositiveInterpolationAlgorithm::InterpolateBoundaries1d:
      interpolation =
          std::make_shared<datadriven::OperationMakePositiveInterpolateBoundaryOfSupport>();
      break;
    case MakePositiveInterpolationAlgorithm::InterpolateFunction:
      interpolation = std::make_shared<datadriven::OperationMakePositiveInterpolateFunction>(f);
      break;
    default:
      interpolation = std::make_shared<datadriven::OperationMakePositiveSetToZero>();
      break;
  }

  // init grid point lists
  addedGridPoints.clear();
  addedGridPointsForPositivity.clear();
}

void OperationMakePositive::makeCurrentNodalValuesPositive(base::Grid& grid,
                                                           base::DataVector& alpha, double tol) {
  // dehierarchize the grid
  auto opHier = op_factory::createOperationHierarchisation(grid);
  opHier->doDehierarchisation(alpha);

  for (size_t i = 0; i < alpha.getSize(); ++i) {
    if (alpha[i] < tol) {
      alpha[i] = 0.0;
    }
  }

  // hierarchize the result
  opHier->doHierarchisation(alpha);
}

void OperationMakePositive::forceNewNodalValuesToBePositive(base::Grid& grid,
                                                            base::DataVector& alpha,
                                                            std::vector<size_t>& newGridPoints,
                                                            double tol) {
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto opEval = op_factory::createOperationEvalNaive(grid);
  base::DataVector x(gridStorage.getDimension());
  for (size_t i : newGridPoints) {
    gridStorage.getCoordinates(gridStorage.getPoint(i), x);
    double nodalValue = opEval->eval(alpha, x);
    if (nodalValue < tol) {
      alpha[i] -= nodalValue;
    }
  }
}

void OperationMakePositive::extractNonExistingCandidatesByLevelSum(
    base::Grid& newGrid, std::vector<std::shared_ptr<base::HashGridPoint>>& candidates,
    size_t currentLevelSum, std::vector<std::shared_ptr<base::HashGridPoint>>& finalCandidates) {
  base::HashGridStorage& gridStorage = newGrid.getStorage();
  for (auto& candidate : candidates) {
    if (!gridStorage.isContaining(*candidate) && candidate->getLevelSum() == currentLevelSum) {
      finalCandidates.push_back(candidate);
    }
  }
}

void OperationMakePositive::addFullGridPoints(
    base::Grid& grid, base::DataVector& alpha,
    std::vector<std::shared_ptr<base::HashGridPoint>>& candidates,
    std::vector<size_t>& addedGridPoints, double tol) {
  base::HashGridStorage& gridStorage = grid.getStorage();
  size_t numDims = gridStorage.getDimension();
  size_t numCandidates = candidates.size();
  base::DataVector x(numDims);

  // compute the function values of the new candidates
  base::DataMatrix data(numCandidates, numDims);
  for (size_t i = 0; i < candidates.size(); ++i) {
    gridStorage.getCoordinates(*candidates[i], x);
    data.setRow(i, x);
  }

  // use streaming eval to support inconsistent grids
  base::DataVector nodalValues(numCandidates);
  std::shared_ptr<datadriven::OperationMultipleEvalConfiguration> evalConfig;
  if (grid.getType() == base::GridType::Linear) {
    evalConfig.reset(new datadriven::OperationMultipleEvalConfiguration(
        datadriven::OperationMultipleEvalType::STREAMING,
        datadriven::OperationMultipleEvalSubType::DEFAULT));
  } else {
    evalConfig.reset(new datadriven::OperationMultipleEvalConfiguration(
        datadriven::OperationMultipleEvalType::DEFAULT));
  }
  std::unique_ptr<base::OperationMultipleEval> opEval(
      op_factory::createOperationMultipleEval(grid, data, *evalConfig));
  opEval->mult(alpha, nodalValues);

  // insert the negative ones to the grid
  size_t oldNumNewGridPointsForPositivity = numAddedGridPointsForPositivity();
  size_t ix = 0;
  for (size_t i = 0; i < candidates.size(); ++i) {
    if (nodalValues[i] < tol) {
      if (generateConsistentGrid) {
        gridStorage.insert(*candidates[i], addedGridPoints);
        ix = gridStorage.getSequenceNumber(*candidates[i]);
      } else {
        ix = gridStorage.insert(*candidates[i]);
        addedGridPoints.push_back(ix);
      }
      addedGridPointsForPositivity.push_back(ix);
    }
  }

  if (verbose) {
    std::cout << "# negative considered : "
              << (numAddedGridPointsForPositivity() - oldNumNewGridPointsForPositivity)
              << " <= " << addedGridPoints.size() << " : # added grid points" << std::endl;
  }
  // add the number of new grid points to the corresponding vector
  size_t numIterations = countAddedGridPointsForPositivityPerIteration.getSize();
  countAddedGridPointsForPositivityPerIteration.resize(numIterations + 1);
  countAddedGridPointsForPositivityPerIteration[numIterations] =
      static_cast<double>(numAddedGridPointsForPositivity() - oldNumNewGridPointsForPositivity);

  // recompute the leaf property
  grid.getStorage().recalcLeafProperty();

  // resize the alpha vector and set the new values to zero
  alpha.resizeZero(gridStorage.getSize());
}

void OperationMakePositive::makePositive(base::Grid& grid, base::DataVector& alpha,
                                         bool forcePositiveNodalValues) {
  // make sure that we use a local basis, otherwise the operation is not supported
  if (grid.getType() != base::GridType::Linear &&
      grid.getType() != base::GridType::LinearBoundary &&
      grid.getType() != base::GridType::LinearL0Boundary &&
      grid.getType() != base::GridType::LinearTruncatedBoundary &&
      grid.getType() != base::GridType::LinearClenshawCurtis &&
      grid.getType() != base::GridType::LinearClenshawCurtisBoundary) {
    throw base::factory_exception(
        "OperationMakePositive::makePositive - this operation not implemented for this grid "
        "type");
  }

  // initialize the operation with the current parameter setting
  initialize(grid, alpha);
  interpolation->initialize(grid);

  if (forcePositiveNodalValues) {
    // force the nodal values of the initial grid to be positive
    makeCurrentNodalValuesPositive(grid, alpha);
  }

  size_t currentLevelSum = minimumLevelSum;
  size_t oldGridSize = grid.getSize();

  std::vector<size_t> currentlyAddedGridPoints;
  std::vector<std::shared_ptr<base::HashGridPoint>> candidates;

  while (currentLevelSum <= maximumLevelSum) {
    if (verbose) {
      std::cout << "current level sum     : " << currentLevelSum << "/" << maximumLevelSum
                << std::endl;
    }

    // search the next candidate set
    candidateSearch->nextCandidates(grid, alpha, currentLevelSum, candidates);

    if (verbose) {
      std::cout << "# candidates found    : " << candidates.size();
    }

    if (candidates.size() > 0) {
      // remove existing candidates from the candidate set and sort the rest by level sum
      std::vector<std::shared_ptr<base::HashGridPoint>> finalCandidates;
      extractNonExistingCandidatesByLevelSum(grid, candidates, currentLevelSum, finalCandidates);

      if (verbose) {
        std::cout << " >= " << finalCandidates.size() << " / " << candidateSearch->numCandidates()
                  << std::endl;
      }

      // add the those candidates for which the sparse grid function is negative
      addFullGridPoints(grid, alpha, finalCandidates, currentlyAddedGridPoints);

      if (verbose) {
        std::cout << "# new grid points     : " << currentlyAddedGridPoints.size() << " -> "
                  << grid.getSize() << std::endl;
      }

      // update the hierarchical coefficients for the new points
      interpolation->computeHierarchicalCoefficients(grid, alpha, currentlyAddedGridPoints);
    } else {
      if (verbose) {
        std::cout << std::endl;
      }
    }

    if (verbose) {
      std::cout << "--------------------------------------------------------" << std::endl;
    }

    // increment the level sum
    ++currentLevelSum;

    // update the global index list of new grid points
    addedGridPoints.insert(addedGridPoints.begin(), currentlyAddedGridPoints.begin(),
                           currentlyAddedGridPoints.end());
    currentlyAddedGridPoints.clear();
  }

  if (verbose) {
    std::cout << "# added grid points   : " << numAddedGridPointsForPositivity() << " + "
              << (numAddedGridPoints() - numAddedGridPointsForPositivity()) << " = "
              << numAddedGridPoints() << " + " << oldGridSize << " = " << grid.getSize()
              << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
  }
}

size_t OperationMakePositive::numAddedGridPoints() { return addedGridPoints.size(); }

size_t OperationMakePositive::numAddedGridPointsForPositivity() {
  return addedGridPointsForPositivity.size();
}

base::DataVector& OperationMakePositive::numAddedGridPointsForPositivityPerIteration() {
  return countAddedGridPointsForPositivityPerIteration;
}

std::vector<size_t>& OperationMakePositive::getAddedGridPoints() { return addedGridPoints; }

std::vector<size_t>& OperationMakePositive::getAddedGridPointsForPositivity() {
  return addedGridPointsForPositivity;
}

OperationMakePositiveCandidateSetAlgorithm& OperationMakePositive::getCandidateSetAlgorithm() {
  return *candidateSearch;
}

} /* namespace datadriven */
} /* namespace sgpp */

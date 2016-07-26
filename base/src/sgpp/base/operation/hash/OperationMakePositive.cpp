// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationMakePositive.hpp>
#include <sgpp/base/operation/hash/OperationMakePositiveCandidateSetAlgorithm.hpp>
#include <sgpp/base/operation/hash/OperationMakePositiveInterpolationAlgorithm.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace base {

OperationMakePositive::OperationMakePositive(
    base::Grid& grid, MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm,
    MakePositiveInterpolationAlgorithm interpolationAlgorithm, bool verbose)
    : grid(grid), minimumLevelSum(0), maximumLevelSum(0), verbose(verbose) {
  // set range of level sums for candidate search
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto numDims = gridStorage.getDimension();
  auto maxLevel = gridStorage.getMaxLevel();
  minimumLevelSum = numDims;
  maximumLevelSum = maxLevel * numDims;

  // set the candidate search algorithm
  switch (candidateSearchAlgorithm) {
    case MakePositiveCandidateSearchAlgorithm::FullGrid:
      candidateSearch = std::make_shared<base::OperationMakePositiveLoadFullGridCandidates>(grid);
      break;
    case MakePositiveCandidateSearchAlgorithm::Intersections:
      candidateSearch =
          std::make_shared<base::OperationMakePositiveFindIntersectionCandidates>(grid);
      break;
    default:
      candidateSearch =
          std::make_shared<base::OperationMakePositiveFindIntersectionCandidates>(grid);
      break;
  }
  candidateSearch->setVerbose(verbose);

  // set the interpolation algorithm
  switch (interpolationAlgorithm) {
    case MakePositiveInterpolationAlgorithm::SetToZero:
      interpolationMethod = std::make_shared<base::OperationMakePositiveSetToZero>();
      break;
    default:
      interpolationMethod = std::make_shared<base::OperationMakePositiveSetToZero>();
      break;
  }
}

OperationMakePositive::~OperationMakePositive() {}

void OperationMakePositive::copyGrid(base::Grid& grid, base::Grid*& newGrid) {
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

void OperationMakePositive::makeCurrentNodalValuesPositive(base::DataVector& alpha, double tol) {
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
                                                            std::vector<size_t>& newGridPoints) {
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto opEval = op_factory::createOperationEval(grid);
  DataVector x(gridStorage.getDimension());
  for (size_t i : newGridPoints) {
    gridStorage.getPoint(i).getStandardCoordinates(x);
    double nodalValue = opEval->eval(alpha, x);
    if (nodalValue < 0.0) {
      alpha[i] -= nodalValue;
    }
  }
}

void OperationMakePositive::extractNonExistingCandidatesByLevelSum(
    std::vector<std::shared_ptr<base::HashGridPoint>>& candidates,
    std::vector<std::shared_ptr<base::HashGridPoint>>& finalCandidates, size_t currentLevelSum) {
  base::HashGridStorage& gridStorage = grid.getStorage();
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
  DataVector x(numDims);

  // check if the function value at the new candidates is negative
  DataMatrix data(numCandidates, numDims);
  for (size_t i = 0; i < candidates.size(); ++i) {
    candidates[i]->getStandardCoordinates(x);
    data.setRow(i, x);
  }

  DataVector nodalValues(numCandidates);
  op_factory::createOperationMultipleEval(grid, data)->mult(alpha, nodalValues);

  // insert the negative ones to the grid
  size_t cntConsidered = 0;
  for (size_t i = 0; i < candidates.size(); ++i) {
    if (nodalValues[i] < tol) {
      ++cntConsidered;
      gridStorage.insert(*candidates[i], addedGridPoints);
    }
  }

  if (verbose) {
    std::cout << "# considered          : " << cntConsidered << " -> " << addedGridPoints.size()
              << " : # new grid points" << std::endl;
  }

  // recompute the leaf property
  grid.getStorage().recalcLeafProperty();

  // resize the alpha vector and set the new values to zero
  alpha.resizeZero(gridStorage.getSize());
}

void OperationMakePositive::makePositive(base::Grid*& newGrid, base::DataVector& newAlpha,
                                         bool resetGrid) {
  // copy the current grid
  if (resetGrid) {
    copyGrid(grid, newGrid);
  }

  // force the nodal values of the initial grid to be positive
  makeCurrentNodalValuesPositive(newAlpha);
  size_t currentLevelSum = minimumLevelSum;

  std::vector<size_t> addedGridPoints;
  std::vector<std::shared_ptr<HashGridPoint>> candidates;

  while (currentLevelSum <= maximumLevelSum) {
    if (verbose) {
      std::cout << "current level sum     : " << currentLevelSum << "/" << maximumLevelSum
                << std::endl;
    }

    // search the next candidate set
    candidateSearch->nextCandidates(*newGrid, newAlpha, currentLevelSum, candidates);

    if (verbose) {
      std::cout << "# candidates found    : " << candidates.size();
    }

    if (candidates.size() > 0) {
      // remove existing candidates from the candidate set and sort the rest by level sum
      std::vector<std::shared_ptr<base::HashGridPoint>> finalCandidates;
      extractNonExistingCandidatesByLevelSum(candidates, finalCandidates, currentLevelSum);

      if (verbose) {
        std::cout << " >= " << finalCandidates.size() << std::endl;
      }

      // add the those candidates for which the sparse grid function is negative
      addFullGridPoints(*newGrid, newAlpha, finalCandidates, addedGridPoints);

      if (verbose) {
        std::cout << "# new grid points     : " << addedGridPoints.size() << " -> "
                  << newGrid->getSize() << std::endl;
      }

      // update the hierarchical coefficients for the new points
      interpolationMethod->computeHierarchicalCoefficients(*newGrid, newAlpha, addedGridPoints);
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
    addedGridPoints.clear();
  }
}

} /* namespace base */
} /* namespace sgpp */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OperationMakePositive.hpp"
#include <sgpp/base/operation/hash/OperationMakePositiveCandidateSetAlgorithm.hpp>
#include <sgpp/base/operation/hash/OperationMakePositiveInterpolationAlgorithm.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <vector>
#include <map>

namespace sgpp {
namespace base {

OperationMakePositive::OperationMakePositive(
    base::Grid& grid, MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm,
    MakePositiveInterpolationAlgorithm interpolationAlgorithm)
    : grid(grid), lastMinimumCandidateLevelSum(-1) {
  // set the candidate search algorithm
  switch (candidateSearchAlgorithm) {
    case MakePositiveCandidateSearchAlgorithm::FullGrid:
      candidateSearch = std::make_shared<base::OperationMakePositiveLoadFullGridCandidates>(grid);
      break;
    case MakePositiveCandidateSearchAlgorithm::Intersections:
      candidateSearch = std::make_shared<base::OperationMakePositiveFindIntersectionCandidates>();
      break;
    default:
      candidateSearch = std::make_shared<base::OperationMakePositiveFindIntersectionCandidates>();
      break;
  }
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

void OperationMakePositive::makePositive(base::Grid*& newGrid, base::DataVector& newAlpha) {
  // copy the current grid
  copyGrid(grid, newGrid);

  // force the nodal values of the initial grid to be positive
  makeCurrentNodalValuesPositive(newAlpha);

  base::HashGridStorage& gridStorage = grid.getStorage();
  auto numDims = gridStorage.getDimension();
  //  auto numGridPoints = gridStorage.getSize();
  auto maxLevel = gridStorage.getMaxLevel();

  size_t minLevelSum = maxLevel + numDims - 1;
  size_t maxLevelSum = maxLevel * numDims;

  std::vector<size_t> addedGridPoints;
  std::vector<HashGridPoint*> candidates;

  while (minLevelSum <= maxLevelSum) {
    // search the next candidate set
    candidateSearch->nextCandidates(*newGrid, newAlpha, minLevelSum, candidates);

    if (candidates.size() > 0) {
      addFullGridPoints(*newGrid, newAlpha, candidates, addedGridPoints);
      interpolationMethod->computeHierarchicalCoefficients(*newGrid, newAlpha, addedGridPoints);
    }

    // increment the level sum
    ++minLevelSum;
  }
}

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

void OperationMakePositive::makeCurrentNodalValuesPositive(base::DataVector& alpha) {
  // dehierarchize the grid
  auto opHier = op_factory::createOperationHierarchisation(grid);
  opHier->doDehierarchisation(alpha);

  std::vector<size_t> negativeGridPoints;
  for (size_t i = 0; i < alpha.getSize(); ++i) {
    if (alpha[i] < 0.0) {
      negativeGridPoints.push_back(i);
      alpha[i] = 0.0;
    }
  }

  if (negativeGridPoints.size() > 0) {
    // hierarchize the result
    opHier->doHierarchisation(alpha);
  }
}

void OperationMakePositive::forceNewNodalValuesToBePositive(base::Grid& grid,
                                                            base::DataVector& alpha,
                                                            std::vector<size_t>& newGridPoints) {
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto opEval = op_factory::createOperationEval(grid);
  DataVector x(gridStorage.getDimension());
  for (size_t i : newGridPoints) {
    gridStorage.getCoordinates(gridStorage.getPoint(i), x);
    double nodalValue = opEval->eval(alpha, x);
    if (nodalValue < 0.0) {
      alpha[i] -= nodalValue;
    }
  }
}

void OperationMakePositive::sortNonExistingCandidatesByLevelSum(
    std::vector<base::HashGridPoint*>& candidates,
    std::map<size_t, std::vector<base::HashGridPoint*> >& sortedCandidates) {
  base::HashGridStorage& gridStorage = grid.getStorage();
  for (auto& candidate : candidates) {
    size_t levelSum = candidate->getLevelSum();
    if (levelSum >= lastMinimumCandidateLevelSum && !gridStorage.isContaining(*candidate)) {
      if (sortedCandidates.find(levelSum) != sortedCandidates.end()) {
        std::vector<base::HashGridPoint*> gridPointvector(1);
        sortedCandidates.insert(std::make_pair(levelSum, gridPointvector));
      }

      sortedCandidates[levelSum].push_back(candidate);
    }
  }
}

void OperationMakePositive::addFullGridPoints(base::Grid& grid, base::DataVector& alpha,
                                              std::vector<base::HashGridPoint*>& candidates,
                                              std::vector<size_t>& addedGridPoints) {
  // remove existing candidates from the candidate set and sort the rest by level sum
  std::map<size_t, std::vector<base::HashGridPoint*> > finalCandidates;
  sortNonExistingCandidatesByLevelSum(candidates, finalCandidates);

  base::HashGridStorage& gridStorage = grid.getStorage();
  auto opEval = op_factory::createOperationEval(grid);
  DataVector x(gridStorage.getDimension());

  while (addedGridPoints.size() == 0) {
    if (finalCandidates.find(lastMinimumCandidateLevelSum) != finalCandidates.end()) {
      auto curCandidates = finalCandidates[lastMinimumCandidateLevelSum];

      for (auto& candidate : curCandidates) {
        gridStorage.getCoordinates(*candidate, x);
        double nodalValue = opEval->eval(alpha, x);
        if (nodalValue < 0.0) {
          gridStorage.insert(*candidate, addedGridPoints);
        }
      }
      ++lastMinimumCandidateLevelSum;
    }

    // recompute the leaf property and return the result
    grid.getStorage().recalcLeafProperty();
  }
}

} /* namespace base */
} /* namespace sgpp */

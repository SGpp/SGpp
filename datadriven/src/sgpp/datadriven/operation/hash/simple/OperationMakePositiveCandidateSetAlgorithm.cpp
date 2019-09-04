// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationMakePositiveCandidateSetAlgorithm.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <set>
#include <vector>
#include <map>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace datadriven {

OperationMakePositiveCandidateSetAlgorithm::OperationMakePositiveCandidateSetAlgorithm(
    size_t maxLevel)
    : iteration(0), maxLevel(maxLevel), verbose(false) {}

OperationMakePositiveCandidateSetAlgorithm::~OperationMakePositiveCandidateSetAlgorithm() {}

void OperationMakePositiveCandidateSetAlgorithm::setVerbose(bool pverbose) { verbose = pverbose; }

void OperationMakePositiveCandidateSetAlgorithm::findNodesWithNegativeCoefficients(
    base::DataVector& alpha, std::vector<size_t>& negativeGridPoints, double tol) {
  for (size_t i = 0; i < alpha.getSize(); ++i) {
    if (alpha[i] < tol) {
      negativeGridPoints.push_back(i);
    }
  }
}

base::DataVector& OperationMakePositiveCandidateSetAlgorithm::numCandidatesPerLevel() {
  return gridPointsPerLevel;
}

size_t OperationMakePositiveCandidateSetAlgorithm::costsComputingCandidates() {
  return static_cast<size_t>(costsPerIteration.sum());
}

base::DataVector&
OperationMakePositiveCandidateSetAlgorithm::costsComputingCandidatesPerIteration() {
  return costsPerIteration;
}

// -------------------------------------------------------------------------------------------
OperationMakePositiveFindIntersectionCandidates::OperationMakePositiveFindIntersectionCandidates(
    size_t maxLevel)
    : OperationMakePositiveCandidateSetAlgorithm(maxLevel) {}

OperationMakePositiveFindIntersectionCandidates::
    ~OperationMakePositiveFindIntersectionCandidates() {}

bool OperationMakePositiveFindIntersectionCandidates::haveOverlappingSupport(
    base::HashGridPoint& gpi, base::HashGridPoint& gpj, size_t dim) {
  size_t leveli = gpi.getLevel(dim), indexi = gpi.getIndex(dim);
  size_t levelj = gpj.getLevel(dim), indexj = gpj.getIndex(dim);

  if (leveli == levelj) return indexi == indexj;

  if (leveli < levelj)
    return gpi.isHierarchicalAncestor(gpj, dim);
  else
    return gpj.isHierarchicalAncestor(gpi, dim);
}

bool OperationMakePositiveFindIntersectionCandidates::haveOverlappingSupport(
    base::HashGridPoint& gpi, base::HashGridPoint& gpj) {
  size_t idim = 0;
  size_t numDims = gpi.getDimension();

  while (idim < numDims && haveOverlappingSupport(gpi, gpj, idim)) ++idim;

  // check whether the supports are overlapping in all dimensions
  return idim == numDims;
}

void OperationMakePositiveFindIntersectionCandidates::computeIntersection(
    base::HashGridPoint& gpi, base::HashGridPoint& gpj, base::HashGridPoint& gpintersection) {
  for (size_t idim = 0; idim < gpi.getDimension(); ++idim) {
    if (gpi.getLevel(idim) > gpj.getLevel(idim)) {
      gpintersection.set(idim, gpi.getLevel(idim), gpi.getIndex(idim));
    } else {
      gpintersection.set(idim, gpj.getLevel(idim), gpj.getIndex(idim));
    }
  }
}

bool OperationMakePositiveFindIntersectionCandidates::compareGridPointsByHash(
    const std::shared_ptr<base::HashGridPoint>& lhs,
    const std::shared_ptr<base::HashGridPoint>& rhs) {
  return lhs->getHash() < rhs->getHash();
}

void OperationMakePositiveFindIntersectionCandidates::initializeCandidates(
    base::Grid& grid, std::vector<size_t>& negativeGridPoints) {
  base::HashGridStorage& gridStorage = grid.getStorage();

  // prepare the intersection list for pair wise interactions
  for (auto& iseq : negativeGridPoints) {
    auto gpi = std::make_shared<base::HashGridPoint>(gridStorage.getPoint(iseq));
    intersections[gpi->getHash()] =
        std::make_shared<std::vector<std::shared_ptr<base::HashGridPoint>>>();
    currentIntersections.insert(gpi);
  }

  // check intersection of two grid points
  for (auto i = currentIntersections.begin(); i != currentIntersections.end(); i++) {
    auto gpi = *i;
    for (auto j = i; j != currentIntersections.end(); j++) {
      auto gpj = *j;
      if (gpi->getHash() != gpj->getHash() && haveOverlappingSupport(*gpi, *gpj) &&
          !gpi->isHierarchicalAncestor(*gpj) && !gpj->isHierarchicalAncestor(*gpi)) {
        intersections[gpi->getHash()]->push_back(gpj);
        intersections[gpj->getHash()]->push_back(gpi);
      }
    }
  }

  // sort the overlapping grid points to speed up intersection search
  for (auto& gps : intersections) {
    std::sort(gps.second->begin(), gps.second->end(), compareGridPointsByHash);
  }

  if (verbose) {
    std::cout << "# intersections (k=1) : " << currentIntersections.size() << std::endl;
  }
}

void OperationMakePositiveFindIntersectionCandidates::findIntersections(
    base::Grid& grid, base::DataVector& alpha, size_t levelSum,
    std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>>& res) {
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto numDims = gridStorage.getDimension();

  // reset the costs vectors
  numCandidatesIteration.setAll(0.0);
  costsPerIteration.setAll(0.0);
  std::unordered_map<size_t, std::shared_ptr<std::vector<std::shared_ptr<base::HashGridPoint>>>>
      intersectionsk;

  // check for intersections of more than two grid points
  for (size_t k = 2; k <= numDims; ++k) {
    if (verbose) {
      std::cout << "# intersections (k=" << k << ") : " << currentIntersections.size();
    }

    size_t minLengthOverlapping = grid.getSize();
    size_t maxLengthOverlapping = 0;
    for (auto& gpi : currentIntersections) {
      auto overlappingGridPoints = intersections[gpi->getHash()];

      for (size_t j = 0; j < overlappingGridPoints->size(); j++) {
        costsPerIteration[k - 2]++;
        auto gpj = (*overlappingGridPoints)[j];

        // find intersection and store it
        auto gpintersection = std::make_shared<base::HashGridPoint>(numDims);
        computeIntersection(*gpi, *gpj, *gpintersection);

        // check if the intersection has already been found
        if ((res.find(gpintersection->getHash()) == res.end()) &&
            !gridStorage.isContaining(*gpintersection)) {
          // store the grid point in the result map
          res[gpintersection->getHash()] = gpintersection;
          numCandidatesIteration[k - 2]++;

          // join the sets for possible intersections searches in the
          // next iteration
          auto iintersections = intersections[gpi->getHash()];
          auto jintersections = intersections[gpj->getHash()];

          // compute intersection of both overlapping candidate list
          std::vector<std::shared_ptr<base::HashGridPoint>> commonIntersections;
          std::set_intersection(iintersections->begin(), iintersections->end(),
                                jintersections->begin(), jintersections->end(),
                                std::back_inserter(commonIntersections), compareGridPointsByHash);

          // store the result if the common candidates overlap with the current intersection
          // initialize the grid points that need to be considered
          auto commonOverlappingIntersections =
              std::make_shared<std::vector<std::shared_ptr<base::HashGridPoint>>>(
                  commonIntersections);
          //          auto commonOverlappingIntersections =
          //              std::make_shared<std::vector<std::shared_ptr<base::HashGridPoint>>>();
          //          for (auto& gpk : commonIntersections) {
          //            if (gpi->getHash() != gpk->getHash() && gpj->getHash() != gpk->getHash() &&
          //                haveOverlappingSupport(*gpk, *gpintersection) &&
          //                !gpk->isHierarchicalAncestor(*gpintersection)) {
          //              commonOverlappingIntersections->push_back(gpk);
          //              intersectionsk[gpk->getHash()] = intersections[gpk->getHash()];
          //            }
          //          }

          // store the list of common overlapping intersections for next iteration
          for (auto& gpk : *commonOverlappingIntersections) {
            intersectionsk[gpk->getHash()] = intersections[gpk->getHash()];
          }

          // store the grid point for the next iteration if there are any overlapping
          // other grid points left
          if (commonOverlappingIntersections->size() > 0) {
            nextIntersections.insert(gpintersection);
            intersectionsk[gpintersection->getHash()] = commonOverlappingIntersections;
            // stats
            minLengthOverlapping =
                std::min(minLengthOverlapping, commonOverlappingIntersections->size());
            maxLengthOverlapping =
                std::max(maxLengthOverlapping, commonOverlappingIntersections->size());
          }
        }
      }
    }

    // clear first intersection queue and swap them
    currentIntersections.clear();
    auto hIntersections = currentIntersections;
    currentIntersections = nextIntersections;
    nextIntersections = hIntersections;
    intersections = intersectionsk;

    if (verbose) {
      if (k < numDims) {
        std::cout << " -> [" << minLengthOverlapping << ", " << maxLengthOverlapping << "]";
      }

      std::cout << ": " << numCandidatesIteration[k - 2] << " -> " << res.size()
                << " (costs = " << costsPerIteration[k - 2] << "/" << costsPerIteration.sum() << ")"
                << std::endl;
    }

    // just stop if there are no more grid points to be considered in the current iteration
    if (currentIntersections.size() == 0) {
      break;
    }
  }
}

void OperationMakePositiveFindIntersectionCandidates::nextCandidates(
    base::Grid& grid, base::DataVector& alpha, size_t levelSum,
    std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) {
  if (iteration == 0) {
    // find all the grid points with negative coefficient
    std::vector<size_t> negativeGridPoints;
    findNodesWithNegativeCoefficients(alpha, negativeGridPoints);

    if (verbose) {
      std::cout << "--------------------------------------------------------" << std::endl;
      std::cout << "# negative candidates : " << negativeGridPoints.size() << "/" << alpha.getSize()
                << std::endl;
    }

    // search for intersections
    size_t numDims = grid.getStorage().getDimension();
    costsPerIteration.resizeZero(numDims - 1);
    numCandidatesIteration.resizeZero(numDims - 1);
    initializeCandidates(grid, negativeGridPoints);
    findIntersections(grid, alpha, levelSum, this->candidates);

    if (verbose) {
      size_t numDims = grid.getStorage().getDimension();
      size_t maxLevel = grid.getStorage().getMaxLevel();
      double numFullGridPoints = std::pow(std::pow(2, maxLevel) - 1, numDims);
      std::cout << "# considered intersect: " << this->candidates.size() << " / "
                << numFullGridPoints << " : full grid points (l = " << maxLevel << ")" << std::endl;
      std::cout << "# comparison costs    : " << costsPerIteration.sum() << std::endl;
      std::cout << "--------------------------------------------------------" << std::endl;
    }

    // update stats: count the grid points for each level sum
    base::GridStorage& gs = grid.getStorage();
    size_t maxLevel = gs.getMaxLevel();
    size_t numLevels = numDims * (maxLevel - 1) + 1;
    gridPointsPerLevel.resize(numLevels);
    gridPointsPerLevel.setAll(0.0);
    for (auto& value : this->candidates) {
      gridPointsPerLevel[value.second->getLevelSum() - numDims]++;
    }
  }

  // increment the iteration
  ++iteration;

  // load the candidates with the desired level sum
  candidates.clear();
  for (auto& value : this->candidates) {
    auto gp = value.second;
    if (gp->getLevelSum() == levelSum) {
      candidates.push_back(gp);
    }
  }
}

size_t OperationMakePositiveFindIntersectionCandidates::numCandidates() {
  return candidates.size();
}

base::DataVector& OperationMakePositiveFindIntersectionCandidates::numCandidatesPerIteration() {
  return numCandidatesIteration;
}

// -------------------------------------------------------------------------------------------

OperationMakePositiveFindIntersectionCandidatesJoin::
    OperationMakePositiveFindIntersectionCandidatesJoin(size_t maxLevel)
    : OperationMakePositiveFindIntersectionCandidates(maxLevel) {}

OperationMakePositiveFindIntersectionCandidatesJoin::
    ~OperationMakePositiveFindIntersectionCandidatesJoin() {}

void OperationMakePositiveFindIntersectionCandidatesJoin::initializeCandidates(
    base::Grid& grid, std::vector<size_t>& negativeGridPoints) {
  base::HashGridStorage& gridStorage = grid.getStorage();

  // prepare the intersection list for pair wise interactions
  for (auto& iseq : negativeGridPoints) {
    auto gpi = std::make_shared<base::HashGridPoint>(gridStorage.getPoint(iseq));
    intersections[gpi->getHash()] =
        std::make_shared<std::vector<std::shared_ptr<base::HashGridPoint>>>();
    currentIntersections.insert(gpi);
  }

  // check intersection of two grid points
  for (auto i = currentIntersections.begin(); i != currentIntersections.end(); i++) {
    auto gpi = *i;
    for (auto j = i; j != currentIntersections.end(); j++) {
      auto gpj = *j;
      if (gpi->getHash() != gpj->getHash() && haveOverlappingSupport(*gpi, *gpj) &&
          !gpi->isHierarchicalAncestor(*gpj) && !gpj->isHierarchicalAncestor(*gpi)) {
        intersections[gpi->getHash()]->push_back(gpj);
      }
    }
  }

  // sort the overlapping grid points to speed up intersection search
  for (auto& gps : intersections) {
    std::sort(gps.second->begin(), gps.second->end(), compareGridPointsByHash);
  }

  if (verbose) {
    std::cout << "# intersections (k=1) : " << currentIntersections.size() << std::endl;
  }
}

void OperationMakePositiveFindIntersectionCandidatesJoin::findIntersections(
    base::Grid& grid, base::DataVector& alpha, size_t levelSum,
    std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>>& res) {
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto numDims = gridStorage.getDimension();

  // reset the costs vectors
  numCandidatesIteration.setAll(0.0);
  costsPerIteration.setAll(0.0);
  std::unordered_map<size_t, std::shared_ptr<std::vector<std::shared_ptr<base::HashGridPoint>>>>
      intersectionsk;

  // check for intersections of more than two grid points
  for (size_t k = 2; k <= numDims; ++k) {
    if (verbose) {
      std::cout << "# intersections (k=" << k << ") : " << currentIntersections.size();
    }

    size_t minLengthOverlapping = grid.getSize();
    size_t maxLengthOverlapping = 0;
    for (auto& gpi : currentIntersections) {
      auto& overlappingGridPoints = intersections[gpi->getHash()];

      for (size_t j = 0; j < overlappingGridPoints->size(); j++) {
        costsPerIteration[k - 2]++;
        auto& gpj = (*overlappingGridPoints)[j];

        // find intersection and store it
        auto gpintersection = std::make_shared<base::HashGridPoint>(numDims);
        computeIntersection(*gpi, *gpj, *gpintersection);

        if (nextIntersections.find(gpintersection) != nextIntersections.end()) {
          continue;
        }

        // check if the intersection has already been found
        if ((res.find(gpintersection->getHash()) == res.end()) &&
            !gridStorage.isContaining(*gpintersection)) {
          // store the grid point in the result map
          res[gpintersection->getHash()] = gpintersection;
          numCandidatesIteration[k - 2]++;
        }

        // store intersection with updated overlapping grid points for next iteration
        auto commonOverlappingIntersections =
            std::make_shared<std::vector<std::shared_ptr<base::HashGridPoint>>>();
        for (auto& gpk : *intersections[gpj->getHash()]) {
          if (haveOverlappingSupport(*gpk, *gpintersection) &&
              !gpk->isHierarchicalAncestor(*gpintersection) &&
              !gpintersection->isHierarchicalAncestor(*gpk)) {
            commonOverlappingIntersections->push_back(gpk);
            intersectionsk[gpk->getHash()] = intersections[gpk->getHash()];
          }
        }

        // store the grid point for the next iteration if there are any overlapping
        // other grid points left
        if (commonOverlappingIntersections->size() > 0) {
          nextIntersections.insert(gpintersection);
          intersectionsk[gpintersection->getHash()] = commonOverlappingIntersections;
          // stats
          minLengthOverlapping =
              std::min(minLengthOverlapping, commonOverlappingIntersections->size());
          maxLengthOverlapping =
              std::max(maxLengthOverlapping, commonOverlappingIntersections->size());
        }
      }
    }

    // clear first intersection queue and swap them
    currentIntersections.clear();
    auto hIntersections = currentIntersections;
    currentIntersections = nextIntersections;
    nextIntersections = hIntersections;
    intersections = intersectionsk;

    if (verbose) {
      if (k < numDims) {
        std::cout << " -> [" << minLengthOverlapping << ", " << maxLengthOverlapping << "]";
      }

      std::cout << ": " << numCandidatesIteration[k - 2] << " -> " << res.size()
                << " (costs = " << costsPerIteration[k - 2] << "/" << costsPerIteration.sum() << ")"
                << std::endl;
    }

    // just stop if there are no more grid points to be considered in the current iteration
    if (currentIntersections.size() == 0) {
      break;
    }
  }
}

// -------------------------------------------------------------------------------------------
OperationMakePositiveLoadFullGridCandidates::OperationMakePositiveLoadFullGridCandidates(
    size_t maxLevel)
    : OperationMakePositiveCandidateSetAlgorithm(maxLevel) {}

OperationMakePositiveLoadFullGridCandidates::~OperationMakePositiveLoadFullGridCandidates() {}

void OperationMakePositiveLoadFullGridCandidates::initializeFullGrid(base::Grid& grid) {
  size_t numDims = grid.getStorage().getDimension();
  fullGrid.reset(base::Grid::createLinearGrid(numDims));
  fullGrid->getGenerator().full(maxLevel);
}

void OperationMakePositiveLoadFullGridCandidates::nextCandidates(
    base::Grid& grid, base::DataVector& alpha, size_t levelSum,
    std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) {
  if (iteration == 0) {
    initializeFullGrid(grid);

    base::HashGridStorage& fullGridStorage = fullGrid->getStorage();
    size_t numDims = fullGridStorage.getDimension();
    size_t maxLevel = fullGridStorage.getMaxLevel();

    // update stats
    candidatesPerIteration.resize(1);
    candidatesPerIteration[0] = static_cast<double>(fullGridStorage.getSize());

    // count the grid points for each level sum
    size_t numLevels = numDims * (maxLevel - 1) + 1;
    gridPointsPerLevel.resize(numLevels);
    gridPointsPerLevel.setAll(0.0);
    for (size_t i = 0; i < fullGridStorage.getSize(); ++i) {
      auto gp = std::make_shared<base::HashGridPoint>(fullGridStorage.getPoint(i));
      gridPointsPerLevel[gp->getLevelSum() - numDims]++;
    }
  }

  base::HashGridStorage& fullGridStorage = fullGrid->getStorage();
  candidates.clear();
  for (size_t i = 0; i < fullGridStorage.getSize(); ++i) {
    auto gp = std::make_shared<base::HashGridPoint>(fullGridStorage.getPoint(i));
    if (gp->getLevelSum() == levelSum) {
      candidates.push_back(gp);
    }
  }

  iteration++;
}

size_t OperationMakePositiveLoadFullGridCandidates::numCandidates() { return fullGrid->getSize(); }

base::DataVector& OperationMakePositiveLoadFullGridCandidates::numCandidatesPerIteration() {
  return candidatesPerIteration;
}

// -------------------------------------------------------------------------------------------

OperationMakePositiveHybridFindIntersectionCandidates::
    OperationMakePositiveHybridFindIntersectionCandidates(size_t maxLevel)
    : OperationMakePositiveFindIntersectionCandidates(maxLevel) {}

OperationMakePositiveHybridFindIntersectionCandidates::
    ~OperationMakePositiveHybridFindIntersectionCandidates() {}

void OperationMakePositiveHybridFindIntersectionCandidates::findIntersections(
    base::Grid& grid, base::DataVector& alpha, size_t levelSum,
    std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>>& res) {
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto numDims = gridStorage.getDimension();
  // resize stats arrays
  costsPerIteration.resizeZero(iteration + 1);
  numCandidatesIteration.resizeZero(iteration + 1);

  // check for intersections of more than two grid points
  std::set<std::shared_ptr<base::HashGridPoint>, HashGridPointCompare> considerInNextIteration;

  std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEvalNaive(grid));
  base::DataVector x(numDims);
  // just stop if there are no more grid points to be considered in the current iteration
  while (currentIntersections.size() > 0) {
    // make sure that no grid points have a larger level sum than the current one
    for (auto& candidate : currentIntersections) {
      if (candidate->getLevelSum() < levelSum) {
        nextIntersections.insert(candidate);
      } else {
        considerInNextIteration.insert(candidate);
      }
    }
    currentIntersections.clear();
    currentIntersections.swap(nextIntersections);

    if (currentIntersections.size() == 0) {
      break;
    }

    size_t minLengthOverlapping = grid.getSize();
    size_t maxLengthOverlapping = 0;
    for (auto& gpi : currentIntersections) {
      auto& overlappingGridPoints = intersections[gpi->getHash()];

      for (size_t j = 0; j < overlappingGridPoints->size(); j++) {
        costsPerIteration[iteration]++;
        auto& gpj = (*overlappingGridPoints)[j];

        // find intersection and store it
        auto gpintersection = std::make_shared<base::HashGridPoint>(numDims);
        computeIntersection(*gpi, *gpj, *gpintersection);

        if (nextIntersections.find(gpintersection) != nextIntersections.end()) {
          continue;
        }

        // check if the intersection has already been found
        if ((res.find(gpintersection->getHash()) == res.end()) &&
            !gridStorage.isContaining(*gpintersection)) {
          // store the grid point in the result map
          res[gpintersection->getHash()] = gpintersection;
          numCandidatesIteration[iteration]++;
        }

        // store intersection with updated overlapping grid points for next iteration
        auto commonOverlappingIntersections =
            std::make_shared<std::vector<std::shared_ptr<base::HashGridPoint>>>();
        for (auto& gpk : *intersections[gpj->getHash()]) {
          if (haveOverlappingSupport(*gpk, *gpintersection) &&
              !gpk->isHierarchicalAncestor(*gpintersection) &&
              !gpintersection->isHierarchicalAncestor(*gpk)) {
            commonOverlappingIntersections->push_back(gpk);
          }
        }

        // store the grid point for the next iteration if there are any overlapping
        // other grid points left
        if (commonOverlappingIntersections->size() > 0) {
          nextIntersections.insert(gpintersection);
          intersections[gpintersection->getHash()] = commonOverlappingIntersections;
          // stats
          minLengthOverlapping =
              std::min(minLengthOverlapping, commonOverlappingIntersections->size());
          maxLengthOverlapping =
              std::max(maxLengthOverlapping, commonOverlappingIntersections->size());
        }
      }
    }

    // clear first intersection queue and swap them
    currentIntersections.clear();
    auto hIntersections = currentIntersections;
    currentIntersections = nextIntersections;
    nextIntersections = hIntersections;

    if (verbose) {
      if (iteration < numDims) {
        std::cout << "iteration = " << iteration << " -> [" << minLengthOverlapping << ", "
                  << maxLengthOverlapping << "]";
      }

      std::cout << ": " << numCandidatesIteration[iteration] << " -> " << res.size()
                << " (costs = " << costsPerIteration[iteration] << "/" << costsPerIteration.sum()
                << ")" << std::endl;
    }
  }

  // copy the remaining intersections
  currentIntersections.insert(considerInNextIteration.begin(), considerInNextIteration.end());
}

void OperationMakePositiveHybridFindIntersectionCandidates::nextCandidates(
    base::Grid& grid, base::DataVector& alpha, size_t levelSum,
    std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) {
  if (verbose) {
    std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
  }

  size_t numDims = grid.getStorage().getDimension();

  if (iteration == 0) {
    // find all full grid points with |l|_\infty <= fullGridLevel
    std::vector<size_t> negativeGridPoints;
    findNodesWithNegativeCoefficients(alpha, negativeGridPoints);

    if (verbose) {
      std::cout << "# negative candidates : " << negativeGridPoints.size() << "/" << alpha.getSize()
                << std::endl;
    }

    // search for intersections
    initializeCandidates(grid, negativeGridPoints);
  }

  // find the necessary intersections for the current level sum
  findIntersections(grid, alpha, levelSum, this->candidates);

  if (verbose) {
    size_t maxLevel = grid.getStorage().getMaxLevel();
    double numFullGridPoints = std::pow(std::pow(2, maxLevel) - 1, numDims);

    std::cout << "# considered intersect: " << this->candidates.size() << " / " << numFullGridPoints
              << " : full grid points (l = " << maxLevel << ")" << std::endl;
    std::cout << "# comparison costs    : " << costsPerIteration.sum() << std::endl;
    std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
  }

  // increment the iteration
  ++iteration;

  // load the candidates with the desired level sum
  candidates.clear();
  for (auto& value : this->candidates) {
    auto gp = value.second;
    if (gp->getLevelSum() == levelSum) {
      candidates.push_back(gp);
    }
  }
}

} /* namespace datadriven */
} /* namespace sgpp */

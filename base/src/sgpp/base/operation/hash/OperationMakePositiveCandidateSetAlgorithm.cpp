// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OperationMakePositiveCandidateSetAlgorithm.hpp"
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <vector>
#include <map>
#include <algorithm>

namespace sgpp {
namespace base {

OperationMakePositiveCandidateSetAlgorithm::OperationMakePositiveCandidateSetAlgorithm() {}

OperationMakePositiveCandidateSetAlgorithm::~OperationMakePositiveCandidateSetAlgorithm() {}

void OperationMakePositiveCandidateSetAlgorithm::findNodesWithNegativeCoefficients(
    base::DataVector& alpha, std::vector<size_t>& negativeGridPoints) {
  for (size_t i = 0; i < alpha.getSize(); ++i) {
    if (alpha[i] < 0.0) {
      negativeGridPoints.push_back(i);
    }
  }
}

// -------------------------------------------------------------------------------------------
OperationMakePositiveFindIntersectionCandidates::OperationMakePositiveFindIntersectionCandidates(
    base::Grid& grid)
    : iteration(0), verbose(true) {}

OperationMakePositiveFindIntersectionCandidates::
    ~OperationMakePositiveFindIntersectionCandidates() {}

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

void OperationMakePositiveFindIntersectionCandidates::initializeCandidates(
    base::Grid& grid, std::vector<size_t>& negativeGridPoints) {
  base::HashGridStorage& gridStorage = grid.getStorage();

  // prepare the intersection list for pair wise interactions
  for (auto& iseq : negativeGridPoints) {
    // set the first element to be the current grid point
    auto gpi = std::make_shared<HashGridPoint>(gridStorage.getPoint(iseq));
    intersections[gpi->getHash()] = std::make_shared<std::vector<std::shared_ptr<HashGridPoint>>>();
    intersections[gpi->getHash()]->push_back(gpi);
    currentIntersections.insert(std::make_pair(gpi->getHash(), gpi));
  }

  // check intersection of two grid points
  for (size_t i = 0; i < negativeGridPoints.size(); ++i) {
    auto iseq = negativeGridPoints[i];
    auto gpi = std::make_shared<HashGridPoint>(gridStorage.getPoint(iseq));
    for (size_t j = i + 1; j < negativeGridPoints.size(); ++j) {
      auto jseq = negativeGridPoints[j];
      auto gpj = std::make_shared<HashGridPoint>(gridStorage.getPoint(jseq));
      if (gpi->hasOverlappingSupport(*gpj) ||
          (!gpi->isHierarchicalAncestor(*gpj) && !gpj->isHierarchicalAncestor(*gpi))) {
        intersections[gpi->getHash()]->push_back(gpj);
        intersections[gpj->getHash()]->push_back(gpi);
      }
    }
  }

  if (verbose) {
    std::cout << "# intersections (k=1) : " << currentIntersections.size() << std::endl;
  }
}

void OperationMakePositiveFindIntersectionCandidates::findIntersections(
    base::Grid& grid, size_t levelSum, std::map<size_t, std::shared_ptr<HashGridPoint>>& res) {
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto numDims = gridStorage.getDimension();

  // check for intersections of more than two grid points
  std::map<size_t, std::shared_ptr<HashGridPoint>> intersectionsForNextIteration;
  for (size_t k = 2; k <= numDims; ++k) {
    for (auto& candidates : currentIntersections) {
      auto overlappingGridPoints = *intersections[candidates.first];
      auto gpi = *overlappingGridPoints[0];

      // check if we need to consider this grid point in the current iteration or not
      if (gpi.getLevelSum() >= levelSum) {
        intersectionsForNextIteration.insert(candidates);
      } else {
        for (size_t j = 1; j < overlappingGridPoints.size(); j++) {
          auto gpj = *overlappingGridPoints[j];

          // find intersection and store it
          auto gpintersection = std::make_shared<HashGridPoint>(numDims);
          computeIntersection(gpi, gpj, *gpintersection);

          // check if the intersection has already been found
          if ((res.find(gpintersection->getHash()) == res.end()) &&
              !gridStorage.isContaining(*gpintersection)) {
            // store the grid point for the next iteration
            auto p = std::make_pair(gpintersection->getHash(), gpintersection);
            res.insert(p);

            // just consider each new intersection once
            if (nextIntersections.find(gpintersection->getHash()) == nextIntersections.end()) {
              nextIntersections.insert(p);
            }

            // join the sets for possible intersections searches in the
            // next iteration
            auto iintersections = intersections[gpi.getHash()];
            auto jintersections = intersections[gpj.getHash()];

            // initialize the grid points that need to be considered
            auto vec = std::make_shared<std::vector<std::shared_ptr<HashGridPoint>>>(1);
            (*vec)[0] = gpintersection;
            intersections[gpintersection->getHash()] = vec;

            // compute intersection of both overlapping candidate list
            std::vector<std::shared_ptr<HashGridPoint>> commonIntersections;
            std::sort(iintersections->begin(), iintersections->end());
            std::sort(jintersections->begin(), jintersections->end());
            std::set_intersection(iintersections->begin(), iintersections->end(),
                                  jintersections->begin(), jintersections->end(),
                                  std::back_inserter(commonIntersections));

            // store the result if the common candidates overlap with the current intersection
            for (auto& gpk : commonIntersections) {
              if (gpk->hasOverlappingSupport(*gpintersection)) {
                intersections[gpintersection->getHash()]->push_back(gpk);
              }
            }
          }
        }
      }
    }

    // clear first intersection queue and swap them
    currentIntersections.clear();
    auto hIntersections = currentIntersections;
    currentIntersections = nextIntersections;
    nextIntersections = hIntersections;

    if (verbose) {
      std::cout << "# intersections (k=" << k << ") : " << currentIntersections.size() << " -> "
                << res.size() << std::endl;
    }

    // just stop if there are no more grid points to be considered in the current iteration
    if (currentIntersections.size() == 0) {
      break;
    }
  }

  // join the maps which we have to consider in the next iteration
  for (auto& candidate : intersectionsForNextIteration) {
    currentIntersections.insert(candidate);
  }
}

void OperationMakePositiveFindIntersectionCandidates::nextCandidates(
    base::Grid& grid, base::DataVector& alpha, size_t levelSum,
    std::vector<std::shared_ptr<HashGridPoint>>& candidates) {
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
    initializeCandidates(grid, negativeGridPoints);
  }

  findIntersections(grid, levelSum, this->candidates);

  if (verbose) {
    std::cout << "# considered intersect: " << this->candidates.size() << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
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
// -------------------------------------------------------------------------------------------
OperationMakePositiveLoadFullGridCandidates::OperationMakePositiveLoadFullGridCandidates(
    base::Grid& grid) {
  base::HashGridStorage& gridStorage = grid.getStorage();

  fullGrid = Grid::createLinearGrid(gridStorage.getDimension());
  fullGrid->getGenerator().full(gridStorage.getMaxLevel());
}

OperationMakePositiveLoadFullGridCandidates::~OperationMakePositiveLoadFullGridCandidates() {}

void OperationMakePositiveLoadFullGridCandidates::nextCandidates(
    base::Grid& grid, base::DataVector& alpha, size_t levelSum,
    std::vector<std::shared_ptr<HashGridPoint>>& candidates) {
  base::HashGridStorage& fullGridStorage = fullGrid->getStorage();
  candidates.clear();
  for (size_t i = 0; i < fullGridStorage.getSize(); ++i) {
    auto gp = std::make_shared<HashGridPoint>(fullGridStorage.getPoint(i));
    if (gp->getLevelSum() == levelSum) {
      candidates.push_back(gp);
    }
  }
}

} /* namespace base */
} /* namespace sgpp */

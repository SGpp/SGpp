// Copyright (sC) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationMakePositiveCandidateSetAlgorithm.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/exception/factory_exception.hpp>

#include <vector>
#include <map>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace datadriven {

OperationMakePositiveCandidateSetAlgorithm::OperationMakePositiveCandidateSetAlgorithm()
    : iteration(0), verbose(false) {}

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

// -------------------------------------------------------------------------------------------
OperationMakePositiveFindIntersectionCandidates::OperationMakePositiveFindIntersectionCandidates()
    : costs(0) {}

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
    currentIntersections.insert(std::make_pair(gpi->getHash(), gpi));
  }

  // check intersection of two grid points
  size_t cntIntersections = 0;
  for (size_t i = 0; i < negativeGridPoints.size(); ++i) {
    auto iseq = negativeGridPoints[i];
    auto gpi = currentIntersections[gridStorage.getPoint(iseq).getHash()];
    for (size_t j = i + 1; j < negativeGridPoints.size(); ++j) {
      auto jseq = negativeGridPoints[j];
      auto gpj = currentIntersections[gridStorage.getPoint(jseq).getHash()];
      if (haveOverlappingSupport(*gpi, *gpj) && !gpi->isHierarchicalAncestor(*gpj) &&
          !gpj->isHierarchicalAncestor(*gpi)) {
        intersections[gpi->getHash()]->push_back(gpj);
        intersections[gpj->getHash()]->push_back(gpi);
        cntIntersections++;
      }
    }
  }

  // sort the overlapping grid points to speed up intersection search
  for (auto& gps : intersections) {
    std::sort(gps.second->begin(), gps.second->end(), compareGridPointsByHash);
  }

  if (verbose) {
    std::cout << "# intersections (k=1) : " << currentIntersections.size() << "/" << costs
              << std::endl;
  }
}

void OperationMakePositiveFindIntersectionCandidates::findIntersections(
    base::Grid& grid, size_t levelSum,
    std::unordered_map<size_t, std::shared_ptr<base::HashGridPoint>>& res) {
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto numDims = gridStorage.getDimension();

  // check for intersections of more than two grid points
  for (size_t k = 2; k <= numDims; ++k) {
    size_t cntNewIntersections = 0;
    size_t currentCosts = 0;
    for (auto& candidates : currentIntersections) {
      auto gpi = candidates.second;
      auto overlappingGridPoints = intersections[candidates.first];
      for (size_t j = 0; j < overlappingGridPoints->size(); j++) {
        currentCosts++;
        auto gpj = (*overlappingGridPoints)[j];

        // find intersection and store it
        auto gpintersection = std::make_shared<base::HashGridPoint>(numDims);
        computeIntersection(*gpi, *gpj, *gpintersection);

        // check if the intersection has already been found
        if ((res.find(gpintersection->getHash()) == res.end()) &&
            !gridStorage.isContaining(*gpintersection)) {
          // store the grid point in the result map
          res[gpintersection->getHash()] = gpintersection;
          ++cntNewIntersections;

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
              std::make_shared<std::vector<std::shared_ptr<base::HashGridPoint>>>();
          for (auto& gpk : commonIntersections) {
            if (haveOverlappingSupport(*gpk, *gpintersection)) {
              commonOverlappingIntersections->push_back(gpk);
            }
          }

          // store the grid point for the next iteration if there are any overlapping
          // other grid points left
          if (commonOverlappingIntersections->size() > 0) {
            nextIntersections[gpintersection->getHash()] = gpintersection;
            intersections[gpintersection->getHash()] = commonOverlappingIntersections;
          }
        }
      }
    }

    // clear first intersection queue and swap them
    currentIntersections.clear();
    auto hIntersections = currentIntersections;
    currentIntersections = nextIntersections;
    nextIntersections = hIntersections;
    costs += currentCosts;

    if (verbose) {
      std::cout << "# intersections (k=" << k << ") : " << cntNewIntersections << " -> "
                << res.size() << " (costs = " << currentCosts << ")" << std::endl;
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
  if (grid.getType() != base::GridType::Linear) {
    throw base::factory_exception(
        "OperationMakePositiveFindIntersectionCandidates is not implemented for this grid type");
  }

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
    costs = 0;
    initializeCandidates(grid, negativeGridPoints);
    findIntersections(grid, levelSum, this->candidates);

    if (verbose) {
      size_t numDims = grid.getStorage().getDimension();
      size_t maxLevel = grid.getStorage().getMaxLevel();
      double numFullGridPoints = std::pow(std::pow(2, maxLevel) - 1, numDims);

      std::cout << "# considered intersect: " << this->candidates.size() << " / "
                << numFullGridPoints << " : full grid points (l = " << maxLevel << ")" << std::endl;
      std::cout << "# comparison costs    : " << costs << std::endl;
      std::cout << "--------------------------------------------------------" << std::endl;
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

// -------------------------------------------------------------------------------------------
OperationMakePositiveLoadFullGridCandidates::OperationMakePositiveLoadFullGridCandidates() {}

OperationMakePositiveLoadFullGridCandidates::~OperationMakePositiveLoadFullGridCandidates() {}

void OperationMakePositiveLoadFullGridCandidates::initializeFullGrid(base::Grid& grid) {
  size_t numDims = grid.getStorage().getDimension();
  size_t maxLevel = grid.getStorage().getMaxLevel();
  fullGrid.reset(base::Grid::createLinearGrid(numDims));
  fullGrid->getGenerator().full(maxLevel);
}

void OperationMakePositiveLoadFullGridCandidates::nextCandidates(
    base::Grid& grid, base::DataVector& alpha, size_t levelSum,
    std::vector<std::shared_ptr<base::HashGridPoint>>& candidates) {
  if (iteration == 0) {
    initializeFullGrid(grid);
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

} /* namespace datadriven */
} /* namespace sgpp */

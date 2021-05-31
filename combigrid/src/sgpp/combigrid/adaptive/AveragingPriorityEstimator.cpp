// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/adaptive/AveragingPriorityEstimator.hpp>

#include <map>
#include <numeric>
#include <utility>
#include <vector>

namespace sgpp {
namespace combigrid {

AveragingPriorityEstimator::AveragingPriorityEstimator(FullGrid::LevelOccupancy levelOccupancy)
    : levelOccupancy(levelOccupancy) {}

double AveragingPriorityEstimator::estimatePriority(
    const LevelVector& levelVector,
    const std::map<LevelVector, double>& deltasOfDownwardNeighbors) const {
  double sumOfNormDividedByNumberOfPoints = 0.;

  for (const std::pair<const std::vector<unsigned int>, double>& mapEntry :
       deltasOfDownwardNeighbors) {
    // get the number of grid points, assuming grids have boundary points

    const index_t numberOfPoints =
        FullGrid::getNumberOfPointsFromLevel(mapEntry.first, levelOccupancy);
    sumOfNormDividedByNumberOfPoints += mapEntry.second / static_cast<double>(numberOfPoints);
  }

  auto priority =
      sumOfNormDividedByNumberOfPoints / static_cast<double>(deltasOfDownwardNeighbors.size());

  return priority;
}

}  // namespace combigrid
}  // namespace sgpp

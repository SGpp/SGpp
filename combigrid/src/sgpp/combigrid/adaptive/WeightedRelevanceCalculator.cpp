// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/adaptive/WeightedRelevanceCalculator.hpp>

#include <algorithm>
#include <numeric>

namespace sgpp {
namespace combigrid {

WeightedRelevanceCalculator::WeightedRelevanceCalculator(
    double weightDeltaInRelationToNumberOfPoints, FullGrid::LevelOccupancy levelOccupancy)
    : weightDeltaInRelationToNumberOfPoints(weightDeltaInRelationToNumberOfPoints),
      levelOccupancy(levelOccupancy) {}

double WeightedRelevanceCalculator::calculate(const LevelVector& levelVector, double delta) const {
  // get the number of grid points, assuming grids have boundary points
  const index_t numberOfPoints = FullGrid::getNumberOfPointsFromLevel(levelVector, levelOccupancy);

  return std::max(
      weightDeltaInRelationToNumberOfPoints * delta,
      (1. - weightDeltaInRelationToNumberOfPoints) / static_cast<double>(numberOfPoints));
}

}  // namespace combigrid
}  // namespace sgpp

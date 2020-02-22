// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/adaptive/WeightedRelevanceCalculator.hpp>

#include <numeric>

namespace sgpp {
namespace combigrid {

WeightedRelevanceCalculator::WeightedRelevanceCalculator(
    double weightDeltaInRelationToNumberOfPoints)
    : weightDeltaInRelationToNumberOfPoints(weightDeltaInRelationToNumberOfPoints) {}

double WeightedRelevanceCalculator::calculate(
    const LevelVector& levelVector, double delta) const {
  auto numPoints = static_cast<index_t>(1)
                   << std::accumulate(levelVector.begin(), levelVector.end(), 0);
  return std::max(weightDeltaInRelationToNumberOfPoints * delta,
                  (1 - weightDeltaInRelationToNumberOfPoints) / static_cast<double>(numPoints));
}

}  // namespace combigrid
}  // namespace sgpp

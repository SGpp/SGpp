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

double AveragingPriorityEstimator::estimatePriority(const LevelVector& levelVector,
    const std::map<LevelVector, double>& deltasOfDownwardNeighbors) const {
  double sumOfNormDividedByNumberOfPoints = 0.;

  for (const std::pair<const std::vector<unsigned int>, double>& mapEntry :
      deltasOfDownwardNeighbors) {
    const index_t sumOfLevelVector = static_cast<index_t>(1) << std::accumulate(
                                          mapEntry.first.begin(), mapEntry.first.end(), 0);
    sumOfNormDividedByNumberOfPoints += mapEntry.second / sumOfLevelVector;
  }

  return sumOfNormDividedByNumberOfPoints / static_cast<double>(deltasOfDownwardNeighbors.size());
}

}  // namespace combigrid
}  // namespace sgpp

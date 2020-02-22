// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/adaptive/AveragingPriorityEstimator.hpp>

#include <map>
#include <numeric>

namespace sgpp {
namespace combigrid {

double AveragingPriorityEstimator::estimatePriority(const LevelVector& levelVector,
    const std::map<LevelVector, double>& deltasOfDownwardNeighbors) const {
  auto normDividedByNumberOfPoints =
      [](double& accumulateResult,
          const std::pair<const std::vector<unsigned int>, double>& mapEntry) {
        auto sumOfLevelVector = static_cast<index_t>(1) << std::accumulate(
                                    mapEntry.first.begin(), mapEntry.first.end(), 0);
        accumulateResult += mapEntry.second / sumOfLevelVector;
        return accumulateResult;
      };
  auto sumOfNormDividedByNumberOfPoints =
      std::accumulate(deltasOfDownwardNeighbors.begin(), deltasOfDownwardNeighbors.end(), 0.,
                      normDividedByNumberOfPoints);
  return sumOfNormDividedByNumberOfPoints / static_cast<double>(deltasOfDownwardNeighbors.size());
}

}  // namespace combigrid
}  // namespace sgpp

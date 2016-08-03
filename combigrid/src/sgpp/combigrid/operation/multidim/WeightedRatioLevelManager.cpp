// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>

#include <algorithm>

namespace sgpp {
namespace combigrid {

WeightedRatioLevelManager::WeightedRatioLevelManager(double w) : LevelManager(), w(w) {}

WeightedRatioLevelManager::WeightedRatioLevelManager(
    std::shared_ptr<AbstractLevelEvaluator> levelEvaluator, double w)
    : LevelManager(levelEvaluator), w(w) {}

WeightedRatioLevelManager::~WeightedRatioLevelManager() {}

double WeightedRatioLevelManager::computePriority(const MultiIndex& level) {
  auto predecessors = getPredecessors(level);

  if (predecessors.empty()) {
    return 1.0;
  }

  double ret = -1.0;
  MultiIndex levelOne(numDimensions, 0.0);
  double levelOneNorm = levelData->get(levelOne)->norm;

  for (auto& predLevel : predecessors) {
    auto predLevelInfo = levelData->get(predLevel);
    double value1 = (1.0 - w) / static_cast<double>(predLevelInfo->numPoints);
    double value2 = w * predLevelInfo->norm / levelOneNorm;
    ret = std::max(value1, value2);
  }

  return ret;
}

} /* namespace combigrid */
} /* namespace sgpp */

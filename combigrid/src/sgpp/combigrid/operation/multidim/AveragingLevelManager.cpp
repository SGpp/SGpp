// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "AveragingLevelManager.hpp"

namespace sgpp {
namespace combigrid {

double AveragingLevelManager::computePriority(const MultiIndex& level) {
  auto predecessors = getPredecessors(level);

  if (predecessors.empty()) {
    return 1.0;
  }

  double sum = 0.0;

  for (auto& predLevel : predecessors) {  // TODO(holzmudd): num points in previous levels
    auto data = levelData->get(predLevel);
    sum += data->norm;
  }

  return sum / static_cast<double>(predecessors.size());
}

AveragingLevelManager::AveragingLevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator)
    : LevelManager(levelEvaluator) {}

AveragingLevelManager::~AveragingLevelManager() {}

} /* namespace combigrid */
} /* namespace sgpp*/

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <algorithm>
#include <sgpp/combigrid/operation/multidim/MaxLevelManager.hpp>

namespace sgpp {
namespace combigrid {

double MaxLevelManager::computePriority(const MultiIndex& level) {
  auto predecessors = getPredecessors(level);

  if (predecessors.empty()) {
    return 0.0;
  }

  double res = 0;
  for (auto& predLevel : predecessors) {
    if (levelData->containsIndex(predLevel)) {
      auto data = levelData->get(predLevel);
      res = std::max(res, data->norm / static_cast<double>(data->maxNewPoints));
    }
  }
  return res;
}

MaxLevelManager::MaxLevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator)
    : LevelManager(levelEvaluator) {}

MaxLevelManager::MaxLevelManager() : LevelManager() {}

MaxLevelManager::~MaxLevelManager() {}

std::shared_ptr<LevelManager> MaxLevelManager::clone() {
  return std::make_shared<MaxLevelManager>(*this);
}

} /* namespace combigrid */
} /* namespace sgpp*/

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>

namespace sgpp {
namespace combigrid {

double AveragingLevelManager::computePriority(const MultiIndex& level) {
  auto predecessors = getPredecessors(level);

  if (predecessors.empty()) {
    return 0.0;
  }

  double sum = 0.0;
  size_t n = 0;
  for (auto& predLevel : predecessors) {
    if (levelData->containsIndex(predLevel)) {
      auto data = levelData->get(predLevel);
      sum += data->norm / static_cast<double>(data->maxNewPoints);
      n += 1;
    }
  }

  if (n == 0) {
    return 0.0;
  } else {
    return sum / static_cast<double>(n);
  }

  //  double maxNorm = 0.0;
  //  for (auto& predLevel : predecessors) {
  //    if (levelData->containsIndex(predLevel)) {
  //      auto data = levelData->get(predLevel);
  //      double dataNorm = data->norm;
  //      if (dataNorm > maxNorm) {
  //        maxNorm = dataNorm;
  //      }
  //    }
  //  }
}

AveragingLevelManager::AveragingLevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator)
    : LevelManager(levelEvaluator) {}

AveragingLevelManager::AveragingLevelManager() : LevelManager() {}

AveragingLevelManager::~AveragingLevelManager() {}

std::shared_ptr<LevelManager> AveragingLevelManager::clone() {
  return std::make_shared<AveragingLevelManager>(*this);
}

} /* namespace combigrid */
} /* namespace sgpp*/

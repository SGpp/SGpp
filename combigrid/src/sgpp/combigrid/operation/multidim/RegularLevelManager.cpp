// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/RegularLevelManager.hpp>

#include <cmath>

namespace sgpp {
namespace combigrid {

double RegularLevelManager::computePriority(const MultiIndex& level) {
  size_t sum = 0;
  for (size_t i = 0; i < level.size(); ++i) {
    sum += level[i];
  }

  return pow(0.5, sum);
}

RegularLevelManager::RegularLevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator)
    : LevelManager(levelEvaluator) {}

RegularLevelManager::RegularLevelManager() {}

std::shared_ptr<LevelManager> RegularLevelManager::clone() {
  return std::make_shared<RegularLevelManager>(*this);
}

RegularLevelManager::~RegularLevelManager() {}

} /* namespace combigrid */
} /* namespace sgpp */

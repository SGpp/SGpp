// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "LevelHelpers.hpp"

#include <vector>
#include <map>
#include <utility>
#include <algorithm>

namespace sgpp {
namespace combigrid {

LevelInfos::LevelInfos() : counterAdaptive(0) {
  infoOnAddedLevels = std::make_shared<
      std::vector<std::shared_ptr<std::map<MultiIndex, std::shared_ptr<LevelInfo>>>>>();
}

LevelInfos::~LevelInfos() {}

void LevelInfos::incrementCounter() {
  infoOnAddedLevels->push_back(
      std::make_shared<std::map<MultiIndex, std::shared_ptr<LevelInfo>>>());
  counterAdaptive++;
}

void LevelInfos::insert(const MultiIndex &level, std::shared_ptr<LevelInfo> levelInfo) {
  std::pair<MultiIndex, std::shared_ptr<LevelInfo>> item(level, levelInfo);
  (*infoOnAddedLevels)[counterAdaptive - 1]->insert(item);
}

void LevelInfos::maxNormPerIteration(sgpp::base::DataVector &maxNorms) {
  maxNorms.resize(infoOnAddedLevels->size());
  size_t i = 0;
  for (auto &istats : *infoOnAddedLevels) {
    double localMax = 0.0;
    for (auto &item : *istats) {
      auto &levelInfo = item.second;
      localMax = std::max(localMax, std::abs(levelInfo->norm));
    }
    maxNorms[i] = localMax;
    ++i;
  }
}

std::shared_ptr<std::vector<std::shared_ptr<std::map<MultiIndex, std::shared_ptr<LevelInfo>>>>>
LevelInfos::getInfos() {
  return infoOnAddedLevels;
}

} /* namespace combigrid */
} /* namespace sgpp*/

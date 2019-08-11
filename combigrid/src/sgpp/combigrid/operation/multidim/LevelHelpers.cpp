// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/LevelHelpers.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

namespace sgpp {
namespace combigrid {

LevelInfos::LevelInfos() : counterAdaptive(0) {
  infoOnAddedLevels = std::make_shared<RefinementInfosPerStep>();
}

LevelInfos::~LevelInfos() {}

void LevelInfos::incrementCounter() {
  if (infoOnAddedLevels->size() == 0 || !(*infoOnAddedLevels)[counterAdaptive - 1].empty()) {
    std::map<MultiIndex, LevelInfo> dict;
    infoOnAddedLevels->push_back(dict);
    counterAdaptive++;
  }
}

void LevelInfos::insert(const MultiIndex &level, LevelInfo &levelInfo) {
  (*infoOnAddedLevels)[counterAdaptive - 1].insert(
      std::pair<MultiIndex, LevelInfo>(level, levelInfo));
}

void LevelInfos::maxNormPerIteration(sgpp::base::DataVector &maxNorms) {
  maxNorms.resize(infoOnAddedLevels->size());
  size_t i = 0;
  for (auto &istats : *infoOnAddedLevels) {
    double localMax = 0.0;
    for (auto &item : istats) {
      auto &levelInfo = item.second;
      auto absNorm = std::fabs(levelInfo.norm);
      if (localMax < absNorm) {
        localMax = absNorm;
      }
    }
    maxNorms[i] = localMax;
    ++i;
  }
}

std::shared_ptr<RefinementInfosPerStep> LevelInfos::getInfos() { return infoOnAddedLevels; }

} /* namespace combigrid */
} /* namespace sgpp*/

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/base/tools/Printer.hpp>

#include <cassert>

namespace sgpp {
namespace combigrid {

std::vector<LevelVector> getLevelVectorsRecursiveLastDim(const LevelVector& maxLevel,
                                                         const LevelVector& minLevel,
                                                         const LevelVector& prefix) {
  assert(prefix.size() == minLevel.size() - 1);
  assert(minLevel.back() <= maxLevel.back());
  std::vector<LevelVector> oneDRange;
  LevelVector current = prefix;
  current.push_back(0);

  for (level_t l = minLevel.back(); l <= maxLevel.back(); ++l) {
    current[prefix.size()] = l;
    oneDRange.push_back(current);
  }

  return oneDRange;
}

std::vector<LevelVector> getLevelVectorsRecursive(const LevelVector& maxLevel,
                                                  const LevelVector& minLevel,
                                                  const LevelVector& prefix) {
  using sgpp::base::operator<<;
  std::vector<LevelVector> currentVector;
  const size_t currentDim = prefix.size();

  // if we are at the last dimension, add one one-dimensional row of the hypercube
  if (currentDim == minLevel.size() - 1) {
    const std::vector<LevelVector> newPole =
        getLevelVectorsRecursiveLastDim(maxLevel, minLevel, prefix);
    currentVector.insert(currentVector.end(), newPole.begin(), newPole.end());
  } else {
    LevelVector newPrefix = prefix;
    newPrefix.push_back(0);

    // else, recurse to the next dimension
    for (level_t l = minLevel[currentDim]; l <= maxLevel[currentDim]; ++l) {
      newPrefix[prefix.size()] = l;
      const std::vector<LevelVector> newHypercube =
          getLevelVectorsRecursive(maxLevel, minLevel, newPrefix);
      currentVector.insert(currentVector.end(), newHypercube.begin(), newHypercube.end());
    }
  }

  return currentVector;
}

std::vector<LevelVector> hyperCubeOfLevelVectors(const LevelVector& maxLevel,
                                                 const LevelVector& minLevel) {
  return getLevelVectorsRecursive(maxLevel, minLevel, LevelVector{});
}

}  // namespace combigrid
}  // namespace sgpp

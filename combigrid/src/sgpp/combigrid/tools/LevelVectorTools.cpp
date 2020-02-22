// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/tools/LevelVectorTools.hpp>

#include <algorithm>
#include <cassert>

namespace sgpp {
namespace combigrid {

std::vector<LevelVector> LevelVectorTools::generateHyperCube(
    const LevelVector& maxLevel, const LevelVector& minLevel) {
  return generateHyperCubeRecursive(maxLevel, minLevel, LevelVector{});
}

std::vector<LevelVector> LevelVectorTools::generateHyperCubeRecursiveLastDim(
    const LevelVector& maxLevel, const LevelVector& minLevel, const LevelVector& prefix) {
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

std::vector<LevelVector> LevelVectorTools::generateHyperCubeRecursive(
    const LevelVector& maxLevel, const LevelVector& minLevel, const LevelVector& prefix) {
  std::vector<LevelVector> currentVector;
  const size_t currentDim = prefix.size();

  // if we are at the last dimension, add one one-dimensional row of the hypercube
  if (currentDim == minLevel.size() - 1) {
    const std::vector<LevelVector> newPole =
        generateHyperCubeRecursiveLastDim(maxLevel, minLevel, prefix);
    currentVector.insert(currentVector.end(), newPole.begin(), newPole.end());
  } else {
    LevelVector newPrefix = prefix;
    newPrefix.push_back(0);

    // else, recurse to the next dimension
    for (level_t l = minLevel[currentDim]; l <= maxLevel[currentDim]; ++l) {
      newPrefix[prefix.size()] = l;
      const std::vector<LevelVector> newHypercube =
          generateHyperCubeRecursive(maxLevel, minLevel, newPrefix);
      currentVector.insert(currentVector.end(), newHypercube.begin(), newHypercube.end());
    }
  }

  return currentVector;
}

std::vector<LevelVector> LevelVectorTools::generateDiagonalWithBoundary(size_t dim, level_t n) {
  if (dim == 0) {
    return std::vector<LevelVector>();
  } else if (dim == 1) {
    return std::vector<LevelVector>{LevelVector{n}};
  } else if (n == 0) {
    return std::vector<LevelVector>{LevelVector(dim, 0)};
  } else {
    std::vector<LevelVector> result;

    for (level_t l = 0; l <= n; l++) {
      for (LevelVector& level : generateDiagonalWithBoundary(dim - 1, n - l)) {
        level.push_back(l);
        result.push_back(level);
      }
    }

    return result;
  }
}

std::vector<LevelVector> LevelVectorTools::generateDiagonalWithoutBoundary(size_t dim, level_t n) {
  if ((dim == 0) || (n < dim)) {
    return std::vector<LevelVector>();
  } else if (dim == 1) {
    return std::vector<LevelVector>{LevelVector{n}};
  } else {
    std::vector<LevelVector> result;

    for (level_t l = 1; l <= n - dim + 1; l++) {
      for (LevelVector& level : generateDiagonalWithoutBoundary(dim - 1, n - l)) {
        level.push_back(l);
        result.push_back(level);
      }
    }

    return result;
  }
}

std::vector<LevelVector> LevelVectorTools::makeDownwardClosed(
    const std::vector<LevelVector>& subspaceLevels, LevelVector lowestLevelVector) {
  assert(lowestLevelVector.size() == subspaceLevels[0].size());
  std::vector<LevelVector> downwardClosedSet = subspaceLevels;

  // for each subspace level, ...
  for (const LevelVector& subspaceLevel : subspaceLevels) {
    // add the full hypercube of lower levels, if not already present
    for (const LevelVector& level : generateHyperCube(subspaceLevel, lowestLevelVector)) {
      if (std::find(downwardClosedSet.begin(), downwardClosedSet.end(), level) ==
          downwardClosedSet.end()) {
        downwardClosedSet.push_back(level);
      }
    }
  }

  std::sort(downwardClosedSet.begin(), downwardClosedSet.end());
  return downwardClosedSet;
}

}  // namespace combigrid
}  // namespace sgpp

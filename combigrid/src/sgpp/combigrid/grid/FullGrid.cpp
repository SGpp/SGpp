// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/FullGrid.hpp>

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

namespace sgpp {
namespace combigrid {

static index_t getNumberOfPointsOnHierarchicalLevel(level_t level,
                                                    FullGrid::LevelOccupancy levelOccupancy) {
  index_t n = 0;
  switch (levelOccupancy) {
    case FullGrid::LevelOccupancy::TwoToThePowerOfL:
      if (level == 0) {
        n = 2;
      } else {
        n = static_cast<index_t>(1) << (level - 1);
      }
    break;
    case FullGrid::LevelOccupancy::Linear:
      throw sgpp::base::not_implemented_exception();
  }
  return n;
}

index_t FullGrid::getNumberOfPointsFromLevel(level_t level, FullGrid::LevelOccupancy levelOccupancy,
                                             bool hasBoundary) {
  auto numberOfPoints = getNumberOfPointsOnHierarchicalLevel(level, levelOccupancy);
  for (level_t l = level - 1; l >= (hasBoundary) ? 0 : 1; ++l) {
    numberOfPoints += getNumberOfPointsOnHierarchicalLevel(l, levelOccupancy);
  }
  return numberOfPoints;
}

index_t FullGrid::getNumberOfPointsFromLevel(const LevelVector& level,
                                             FullGrid::LevelOccupancy levelOccupancy,
                                             bool hasBoundary) {
  std::vector<index_t> pointsInEachDimension;
  std::transform(
      level.begin(), level.end(), std::back_inserter(pointsInEachDimension),
      [levelOccupancy](index_t l) { return getNumberOfPointsFromLevel(l, levelOccupancy); });
  return std::accumulate(pointsInEachDimension.begin(), pointsInEachDimension.end(), 1,
                         std::multiplies<index_t>());
}

bool FullGrid::findGridPointInFullGrid(const base::GridPoint& gridPoint, IndexVector& index) const {
  const LevelVector& levelFullGrid = this->getLevel();
  bool isContained = true;

  for (size_t d = 0; d < gridPoint.getDimension(); d++) {
    const level_t curLevel = gridPoint.getLevel(d);

    if ((curLevel <= levelFullGrid[d]) && (this->hasBoundary() || (curLevel >= 1))) {
      index[d] = gridPoint.getIndex(d) << (levelFullGrid[d] - curLevel);
    } else {
      isContained = false;
      break;
    }
  }

  return isContained;
}

}  // namespace combigrid
}  // namespace sgpp

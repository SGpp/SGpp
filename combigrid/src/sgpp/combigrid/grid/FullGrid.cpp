// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/FullGrid.hpp>

#include <algorithm>
#include <numeric>

namespace sgpp {
namespace combigrid {

index_t FullGrid::getNumberOfPointsOnLevel(level_t level, FullGrid::LevelOccupancy levelOccupancy) {
  switch (levelOccupancy) {
    case FullGrid::LevelOccupancy::TwoToThePowerOfL:
      if (level == 0) {
        return 2;
      } else {
        return static_cast<index_t>(1) << level;
      }
      break;
    case FullGrid::LevelOccupancy::Linear:
    default:
      throw sgpp::base::not_implemented_exception();
      break;
  }
}

index_t FullGrid::getNumberOfPointsOnLevel(const LevelVector& level,
                                           FullGrid::LevelOccupancy levelOccupancy) {
  std::vector<index_t> pointsInEachDimension;
  std::transform(
      level.begin(), level.end(), std::back_inserter(pointsInEachDimension),
      [levelOccupancy](index_t l) { return getNumberOfPointsOnLevel(l, levelOccupancy); });
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

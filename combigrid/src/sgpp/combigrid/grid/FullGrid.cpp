// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/FullGrid.hpp>

namespace sgpp {
namespace combigrid {

bool FullGrid::findGridPointInFullGrid(const base::GridPoint& gridPoint,
    IndexVector& index) const {
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
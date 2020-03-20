// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/adaptive/PriorityEstimator.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>

#include <map>

namespace sgpp {
namespace combigrid {

/**
 * @brief the AveragingLevelManager from holzmudd's combigrid module:
 * The priority returned is the average of delta divided by the number of grid points (as a measure
 * for work intensity) over all the predecessor grids.
 */
class AveragingPriorityEstimator : public PriorityEstimator {
 public:
  AveragingPriorityEstimator(
      FullGrid::LevelOccupancy levelOccupancy = FullGrid::LevelOccupancy::TwoToThePowerOfL);

  double estimatePriority(
      const LevelVector& levelVector,
      const std::map<LevelVector, double>& deltasOfDownwardNeighbors) const override;

 private:
  FullGrid::LevelOccupancy levelOccupancy;
};

}  // namespace combigrid
}  // namespace sgpp

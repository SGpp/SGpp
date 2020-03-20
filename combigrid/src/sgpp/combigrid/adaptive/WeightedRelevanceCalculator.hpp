// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/adaptive/RelevanceCalculator.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>

namespace sgpp {
namespace combigrid {

/**
 * @brief the weighted relevance calculator introduced in [0]
 *
 * [0] Gerstner, T. and Griebel, M., 2003. Dimension–adaptive tensor–product quadrature.
 * Computing, 71(1), pp.65-87.
 */
class WeightedRelevanceCalculator : public RelevanceCalculator {
 public:
  explicit WeightedRelevanceCalculator(
      double weightDeltaInRelationToNumberOfPoints = 1.,
      FullGrid::LevelOccupancy levelOccupancy = FullGrid::LevelOccupancy::TwoToThePowerOfL);

  double calculate(const LevelVector& levelVector, double delta) const override;

 private:
  double weightDeltaInRelationToNumberOfPoints;
  FullGrid::LevelOccupancy levelOccupancy;
};

}  // namespace combigrid
}  // namespace sgpp

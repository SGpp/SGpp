// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/adaptive/RelevanceCalculator.hpp>

namespace sgpp {
namespace combigrid {

/**
 * @brief the weighted relevance calculator introduced in [0]
 */
class WeightedRelevanceCalculator : public RelevanceCalculator {
 public:
  explicit WeightedRelevanceCalculator(double weightDeltaInRelationToNumberOfPoints = 1.);

  double calculate(const LevelVector& levelVector, double delta) const override;

 private:
  double weightDeltaInRelationToNumberOfPoints;
};

}  // namespace combigrid
}  // namespace sgpp

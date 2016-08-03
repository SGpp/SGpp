// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/operation/multidim/LevelManager.hpp>

namespace sgpp {
namespace combigrid {

/**
 * This refinement criterion implements the one proposed by ´Gerstner and Griebel:
 * Dimension-adaptive Tensor-Product Quadrature, p. 9´
 */
class WeightedRatioLevelManager : public LevelManager {
 protected:
  virtual double computePriority(MultiIndex const &level);

 public:
  explicit WeightedRatioLevelManager(double w = 0.8);
  explicit WeightedRatioLevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator,
                                     double w = 0.8);
  virtual ~WeightedRatioLevelManager();

 private:
  double w;
};

} /* namespace combigrid */
} /* namespace sgpp */

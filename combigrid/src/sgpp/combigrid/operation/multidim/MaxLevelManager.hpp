// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/operation/multidim/LevelManager.hpp>

namespace sgpp {
namespace combigrid {

/**
 * This is a simple LevelManager implementation that does level norm prediction for adaptive
 * refinement by using the maximum norm of all predecessor levels divided by their number of grid
 * points
 */
class MaxLevelManager : public LevelManager {
 protected:
  virtual double computePriority(MultiIndex const &level);

 public:
  explicit MaxLevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator);
  MaxLevelManager();

  virtual std::shared_ptr<LevelManager> clone();

  virtual ~MaxLevelManager();
};

} /* namespace combigrid */
} /* namespace sgpp*/

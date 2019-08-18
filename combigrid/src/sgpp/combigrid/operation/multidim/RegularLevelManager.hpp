// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/operation/multidim/LevelManager.hpp>

namespace sgpp {
namespace combigrid {
/**
 * This class provides an adaption strategy that is equivalent to adding levels regularly. It can be
 * used to polymorphically switch between regular and adaptive level generation by using this
 * LevelManager or another LevelManager.
 */
class RegularLevelManager : public LevelManager {
  virtual double computePriority(MultiIndex const &level);

 public:
  explicit RegularLevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator);
  RegularLevelManager();

  virtual std::shared_ptr<LevelManager> clone();

  virtual ~RegularLevelManager();
};

} /* namespace combigrid */
} /* namespace sgpp */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include "LevelManager.hpp"

namespace sgpp {
namespace combigrid {

class RegularLevelManager : public LevelManager {
  virtual double computePriority(MultiIndex const &level);

 public:
  explicit RegularLevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator);
  RegularLevelManager();

  virtual ~RegularLevelManager();
};

} /* namespace combigrid */
} /* namespace sgpp */

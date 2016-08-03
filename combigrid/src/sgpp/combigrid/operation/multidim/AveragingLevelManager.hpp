// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_AVERAGINGLEVELMANAGER_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_AVERAGINGLEVELMANAGER_HPP_

#include "LevelManager.hpp"

namespace sgpp {
namespace combigrid {

class AveragingLevelManager : public LevelManager {
 protected:
  virtual double computePriority(MultiIndex const &level);

 public:
  explicit AveragingLevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator);
  AveragingLevelManager();

  virtual ~AveragingLevelManager();
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_AVERAGINGLEVELMANAGER_HPP_ */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_AVERAGINGLEVELMANAGER_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_AVERAGINGLEVELMANAGER_HPP_

#include <sgpp/combigrid/operation/multidim/LevelManager.hpp>

namespace sgpp {
namespace combigrid {

/**
 * This is a simple LevelManager implementation that does level norm prediction for adaptive
 * refinement by averaging the norms of predecessor levels divided by their number of new points.
 */
class AveragingLevelManager : public LevelManager {
 protected:
  virtual double computePriority(MultiIndex const &level);

 public:
  explicit AveragingLevelManager(std::shared_ptr<AbstractLevelEvaluator> levelEvaluator);
  AveragingLevelManager();

  virtual std::shared_ptr<LevelManager> clone();

  virtual ~AveragingLevelManager();
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_AVERAGINGLEVELMANAGER_HPP_ */

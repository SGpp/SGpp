// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_CUSTOMGROWTHSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_CUSTOMGROWTHSTRATEGY_HPP_

#include <sgpp/combigrid/grid/growth/AbstractGrowthStrategy.hpp>

#include <functional>

namespace sgpp {
namespace combigrid {

/**
 * The level-numPoints-mapping can be customized with a function.
 */
class CustomGrowthStrategy : public AbstractGrowthStrategy {
  std::function<size_t(size_t)> func;

 public:
  explicit CustomGrowthStrategy(std::function<size_t(size_t)> func);
  virtual ~CustomGrowthStrategy();

  virtual size_t numPoints(size_t level);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_CUSTOMGROWTHSTRATEGY_HPP_ */

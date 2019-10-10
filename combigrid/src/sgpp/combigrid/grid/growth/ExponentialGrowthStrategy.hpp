// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_EXPONENTIALGROWTHSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_EXPONENTIALGROWTHSTRATEGY_HPP_

#include <sgpp/combigrid/grid/growth/AbstractGrowthStrategy.hpp>

namespace sgpp {
namespace combigrid {

/**
 * numPoints = (l == 0) ? 1 : 2^l + 1
 * Yields nested points for e. g. uniform points (with boundary) and Clenshaw-Curtis-Points
 */
class ExponentialGrowthStrategy : public AbstractGrowthStrategy {
 public:
  virtual ~ExponentialGrowthStrategy();

  virtual size_t numPoints(size_t level);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_EXPONENTIALGROWTHSTRATEGY_HPP_ */

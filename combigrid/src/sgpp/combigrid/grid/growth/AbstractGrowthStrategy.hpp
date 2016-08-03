// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ABSTRACTGROWTHSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ABSTRACTGROWTHSTRATEGY_HPP_

#include <sgpp/globaldef.hpp>

#include <cstddef>

namespace sgpp {
namespace combigrid {

/**
 * Defines a converter from a level to a number of points, i. e. an abstract base class for
 * level-numPoints mappings.
 * AbstractGrowthStrategy-Objects are used in some subclasses of AbstractPointOrdering.
 */
class AbstractGrowthStrategy {
 public:
  virtual ~AbstractGrowthStrategy();

  virtual size_t numPoints(size_t level) = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ABSTRACTGROWTHSTRATEGY_HPP_ */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ABSTRACTPOINTDISTRIBUTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ABSTRACTPOINTDISTRIBUTION_HPP_

#include <sgpp/globaldef.hpp>

#include <cstddef>

namespace sgpp {
namespace combigrid {

/**
 * An abstract point distribution provides a set of n one-dimensional grid points for each value of
 * n
 */
class AbstractPointDistribution {
 public:
  virtual ~AbstractPointDistribution();

  /**
   * computes the j-th grid point from the set of numPoints grid points.
   */
  virtual double compute(size_t numPoints, size_t j) = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ABSTRACTPOINTDISTRIBUTION_HPP_ */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_UNIFORMPOINTDISTRIBUTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_UNIFORMPOINTDISTRIBUTION_HPP_

#include <sgpp/combigrid/grid/distribution/AbstractPointDistribution.hpp>

namespace sgpp {
namespace combigrid {

/**
 * provides uniform points, i. e. {k/(n-1) for k = 0, ..., n-1} if n >= 2 or {0.5} if n = 1.
 */
class UniformBoundaryPointDistribution : public AbstractPointDistribution {
 public:
  virtual ~UniformBoundaryPointDistribution();

  virtual double compute(size_t numPoints, size_t j);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_UNIFORMPOINTDISTRIBUTION_HPP_ */

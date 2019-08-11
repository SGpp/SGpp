// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_CLENSHAWCURTISDISTRIBUTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_CLENSHAWCURTISDISTRIBUTION_HPP_

#include <sgpp/combigrid/grid/distribution/AbstractPointDistribution.hpp>

namespace sgpp {
namespace combigrid {

/**
 * Provides Clenshaw-Curtis grid points, that is, the points 0.5 - 0.5
 * cos(k*pi/(n-1)) for k = 0, ..., n or 0.5, if n = 1
 */
class ClenshawCurtisDistribution : public AbstractPointDistribution {
 public:
  virtual ~ClenshawCurtisDistribution();

  virtual double compute(size_t numPoints, size_t j);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_CLENSHAWCURTISDISTRIBUTION_HPP_ */

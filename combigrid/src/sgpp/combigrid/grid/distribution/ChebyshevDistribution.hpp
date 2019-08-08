// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SRC_SGPP_COMBIGRID_GRID_DISTRIBUTION_CHEBYSHEVDISTRIBUTION_HPP_
#define SRC_SGPP_COMBIGRID_GRID_DISTRIBUTION_CHEBYSHEVDISTRIBUTION_HPP_

#include <sgpp/combigrid/grid/distribution/AbstractPointDistribution.hpp>

namespace sgpp {
namespace combigrid {

/**
 * Provides Chebyshev grid points, that is, the points
 * 0.5 - 0.5*cos((2k + 1)*pi/(2n)) for k = 0, ..., n-1
 */
class ChebyshevDistribution : public AbstractPointDistribution {
 public:
  virtual ~ChebyshevDistribution();

  virtual double compute(size_t numPoints, size_t j);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* SRC_SGPP_COMBIGRID_GRID_DISTRIBUTION_CHEBYSHEVDISTRIBUTION_HPP_ */

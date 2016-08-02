/*
 * ClenshawCurtisDistribution.hpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_CLENSHAWCURTISDISTRIBUTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_CLENSHAWCURTISDISTRIBUTION_HPP_

#include "AbstractPointDistribution.hpp"

namespace sgpp{
namespace combigrid {

/**
 * Provides Clenshaw-Curtis grid points (Chebyshev grid points), that is, the points 0.5 - 0.5 cos(k*pi/(n-1)) for k = 0, ..., n or 0.5, if n = 1
 */
class ClenshawCurtisDistribution : public AbstractPointDistribution {
public:
	virtual ~ClenshawCurtisDistribution();

	virtual double compute(size_t numPoints, size_t j);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_CLENSHAWCURTISDISTRIBUTION_HPP_ */

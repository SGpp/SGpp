/*
 * UniformPointDistribution.hpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_UNIFORMPOINTDISTRIBUTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_UNIFORMPOINTDISTRIBUTION_HPP_

#include "AbstractPointDistribution.hpp"

namespace SGPP {
namespace combigrid {

/**
 * provides uniform points, i. e. {k/(n-1) for k = 0, ..., n-1} if n >= 2 or {0.5} if n = 1.
 */
class UniformPointDistribution : public AbstractPointDistribution {
public:
	virtual ~UniformPointDistribution();

	virtual SGPP::float_t compute(size_t numPoints, size_t j);
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_UNIFORMPOINTDISTRIBUTION_HPP_ */

/*
 * AbstractPointDistribution.hpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ABSTRACTPOINTDISTRIBUTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ABSTRACTPOINTDISTRIBUTION_HPP_

#include <sgpp/globaldef.hpp>
#include <cstddef>

namespace SGPP {
namespace combigrid {

/**
 * An abstract point distribution provides a set of n one-dimensional grid points for each value of n
 */
class AbstractPointDistribution {
public:
	virtual ~AbstractPointDistribution();

	/**
	 * computes the j-th grid point from the set of numPoints grid points.
	 */
	virtual SGPP::float_t compute(size_t numPoints, size_t j) = 0;
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ABSTRACTPOINTDISTRIBUTION_HPP_ */

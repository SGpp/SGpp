/*
 * ExponentialGrowthStrategy.hpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_EXPONENTIALGROWTHSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_EXPONENTIALGROWTHSTRATEGY_HPP_

#include "AbstractGrowthStrategy.hpp"

namespace sgpp{
namespace combigrid {

/**
 * numPoints = (l == 0) ? 1 : 2^l + 1
 */
class ExponentialGrowthStrategy : public AbstractGrowthStrategy {
public:
	virtual ~ExponentialGrowthStrategy();

	virtual size_t numPoints(size_t level);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_EXPONENTIALGROWTHSTRATEGY_HPP_ */

/*
 * LinearGrowthStrategy.hpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_LINEARGROWTHSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_LINEARGROWTHSTRATEGY_HPP_

#include "AbstractGrowthStrategy.hpp"

namespace sgpp{
namespace combigrid {

/**
 * Calculation: numPoints := 1 + factor * level, where factor is a parameter passed to the class
 */
class LinearGrowthStrategy: public AbstractGrowthStrategy {
	size_t factor;

public:
	LinearGrowthStrategy(size_t factor);
	virtual ~LinearGrowthStrategy();

	virtual size_t numPoints(size_t level);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_LINEARGROWTHSTRATEGY_HPP_ */

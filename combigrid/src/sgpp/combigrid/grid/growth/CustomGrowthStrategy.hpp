/*
 * CustomGrowthStrategy.hpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_CUSTOMGROWTHSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_CUSTOMGROWTHSTRATEGY_HPP_

#include "AbstractGrowthStrategy.hpp"

#include <functional>

namespace sgpp{
namespace combigrid {

/**
 * The level-numPoints-mapping can be customized with a function.
 */
class CustomGrowthStrategy: public AbstractGrowthStrategy {
	std::function<size_t(size_t)> func;

public:
	CustomGrowthStrategy(std::function<size_t(size_t)> func);
	virtual ~CustomGrowthStrategy();

	virtual size_t numPoints(size_t level);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_GROWTH_CUSTOMGROWTHSTRATEGY_HPP_ */

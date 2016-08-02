/*
 * CustomGrowthStrategy.cpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#include "CustomGrowthStrategy.hpp"

namespace sgpp{
namespace combigrid {

CustomGrowthStrategy::CustomGrowthStrategy(std::function<size_t(size_t)> func)
	: func(func) {
}

CustomGrowthStrategy::~CustomGrowthStrategy() {
	// TODO Auto-generated destructor stub
}

size_t CustomGrowthStrategy::numPoints(size_t level) {
	return func(level);
}

} /* namespace combigrid */
} /* namespace sgpp*/

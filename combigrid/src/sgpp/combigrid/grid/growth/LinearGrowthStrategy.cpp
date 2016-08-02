/*
 * LinearGrowthStrategy.cpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#include "LinearGrowthStrategy.hpp"

namespace SGPP {
namespace combigrid {

LinearGrowthStrategy::LinearGrowthStrategy(size_t factor)
	: factor(factor) {
}

LinearGrowthStrategy::~LinearGrowthStrategy() {
}

size_t LinearGrowthStrategy::numPoints(size_t level) {
	return 1 + factor * level;
}

} /* namespace combigrid */
} /* namespace SGPP */

/*
 * ExponentialGrowthStrategy.cpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#include "ExponentialGrowthStrategy.hpp"

namespace SGPP {
namespace combigrid {

ExponentialGrowthStrategy::~ExponentialGrowthStrategy() {
	// TODO Auto-generated destructor stub
}

size_t ExponentialGrowthStrategy::numPoints(size_t level) {
	return level == 0 ? 1 : (static_cast<size_t>(1) << level) + 1;
}

} /* namespace combigrid */
} /* namespace SGPP */

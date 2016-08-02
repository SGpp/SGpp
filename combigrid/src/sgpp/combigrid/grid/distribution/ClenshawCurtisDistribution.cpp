/*
 * ClenshawCurtisDistribution.cpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#include "ClenshawCurtisDistribution.hpp"

#include <cmath>
#include <stdexcept>

namespace sgpp{
namespace combigrid {

ClenshawCurtisDistribution::~ClenshawCurtisDistribution() {
}

double ClenshawCurtisDistribution::compute(size_t numPoints, size_t j) {
	if (j >= numPoints) {
		throw std::logic_error("ClenshawCurtisDistribution::compute: j >= numPoints");
	}

	if (numPoints == 1) {
		return 0.5;
	}

	return 0.5 * (1 - cos((static_cast<double>(j) / static_cast<double>(numPoints - 1)) * M_PI));
}

} /* namespace combigrid */
} /* namespace sgpp*/

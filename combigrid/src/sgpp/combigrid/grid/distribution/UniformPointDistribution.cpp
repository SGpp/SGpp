/*
 * UniformPointDistribution.cpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#include "UniformPointDistribution.hpp"

#include <stdexcept>

namespace SGPP {
namespace combigrid {

UniformPointDistribution::~UniformPointDistribution() {
}

SGPP::float_t UniformPointDistribution::compute(size_t numPoints, size_t j) {
	if(j >= numPoints) {
		throw std::logic_error("UniformPointDistribution::compute: j >= numPoints");
	}

	if(numPoints == 1) {
		return 0.5;
	}

	return static_cast<float_t>(j) / static_cast<float_t>(numPoints-1);
}

} /* namespace combigrid */
} /* namespace SGPP */

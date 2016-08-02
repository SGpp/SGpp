/*
 * AdaptiveRefinementStrategy.cpp
 *
 *  Created on: 08.02.2016
 *      Author: david
 */

#include "AdaptiveRefinementStrategy.hpp"

#include <stdexcept>
#include <cmath>

namespace SGPP {
namespace combigrid {

AdaptiveRefinementStrategy::AdaptiveRefinementStrategy(priority_function func) :
		func(func) {
}

float_t AdaptiveRefinementStrategy::computePriority(std::vector<float_t> const &predecessorNorms, size_t numNewPoints) {
	return func(predecessorNorms, numNewPoints);
}

AdaptiveRefinementStrategy AdaptiveRefinementStrategy::maxStrategy() {
	return AdaptiveRefinementStrategy([](std::vector<float_t> const &predecessorNorms, size_t numNewPoints) -> float_t {
		if(predecessorNorms.size() == 0) throw std::runtime_error("No predecessors in AdaptiveRefinementStrategy");

		float_t val = predecessorNorms[0];

		for(size_t i = 1; i < predecessorNorms.size(); ++i) {
			float_t other = predecessorNorms[i];
			if(other > val) {
				val = other;
			}
		}

		return val / static_cast<float_t>(numNewPoints);
	});
}

AdaptiveRefinementStrategy AdaptiveRefinementStrategy::minStrategy() {
	return AdaptiveRefinementStrategy([](std::vector<float_t> const &predecessorNorms, size_t numNewPoints) -> float_t {
		if(predecessorNorms.size() == 0) throw std::runtime_error("No predecessors in AdaptiveRefinementStrategy");

		float_t val = predecessorNorms[0];

		for(size_t i = 1; i < predecessorNorms.size(); ++i) {
			float_t other = predecessorNorms[i];
			if(other < val) {
				val = other;
			}
		}

		return val / static_cast<float_t>(numNewPoints);
	});
}

AdaptiveRefinementStrategy AdaptiveRefinementStrategy::arithmeticMeanStrategy() {
	return AdaptiveRefinementStrategy([](std::vector<float_t> const &predecessorNorms, size_t numNewPoints) -> float_t {
		if(predecessorNorms.size() == 0) throw std::runtime_error("No predecessors in AdaptiveRefinementStrategy");

		float_t sum = predecessorNorms[0];

		for(size_t i = 1; i < predecessorNorms.size(); ++i) {
			sum += predecessorNorms[i];
		}

		return sum / static_cast<float_t>(numNewPoints);
	});
}

AdaptiveRefinementStrategy AdaptiveRefinementStrategy::geometricMeanStrategy() {
	return AdaptiveRefinementStrategy([](std::vector<float_t> const &predecessorNorms, size_t numNewPoints) -> float_t {
		if(predecessorNorms.size() == 0) throw std::runtime_error("No predecessors in AdaptiveRefinementStrategy");

		float_t prod = predecessorNorms[0];

		for(size_t i = 1; i < predecessorNorms.size(); ++i) {
			prod *= predecessorNorms[i];
		}

		return pow(prod, 1.0 / static_cast<float_t>(numNewPoints));
	});
}

} /* namespace combigrid */
} /* namespace SGPP */

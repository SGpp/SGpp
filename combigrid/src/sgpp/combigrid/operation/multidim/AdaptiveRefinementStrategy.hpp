/*
 * AdaptiveRefinementStrategy.hpp
 *
 *  Created on: 08.02.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ADAPTIVEREFINEMENTSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ADAPTIVEREFINEMENTSTRATEGY_HPP_

#include <sgpp/globaldef.hpp>

#include <vector>
#include <functional>

namespace SGPP {
namespace combigrid {

class AdaptiveRefinementStrategy {
public:
	typedef std::function<float_t(std::vector<float_t> const &, size_t)> priority_function;

	AdaptiveRefinementStrategy(priority_function func);

	float_t computePriority(std::vector<float_t> const &predecessorNorms, size_t numNewPoints);

	static AdaptiveRefinementStrategy maxStrategy();
	static AdaptiveRefinementStrategy minStrategy();
	static AdaptiveRefinementStrategy arithmeticMeanStrategy();
	static AdaptiveRefinementStrategy geometricMeanStrategy();

private:
	priority_function func;
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ADAPTIVEREFINEMENTSTRATEGY_HPP_ */

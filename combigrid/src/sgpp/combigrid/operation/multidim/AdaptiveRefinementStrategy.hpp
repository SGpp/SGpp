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

namespace sgpp{
namespace combigrid {

class AdaptiveRefinementStrategy {
public:
	typedef std::function<double(std::vector<double> const &, size_t)> priority_function;

	AdaptiveRefinementStrategy(priority_function func);

	double computePriority(std::vector<double> const &predecessorNorms, size_t numNewPoints);

	static AdaptiveRefinementStrategy maxStrategy();
	static AdaptiveRefinementStrategy minStrategy();
	static AdaptiveRefinementStrategy arithmeticMeanStrategy();
	static AdaptiveRefinementStrategy geometricMeanStrategy();

private:
	priority_function func;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ADAPTIVEREFINEMENTSTRATEGY_HPP_ */

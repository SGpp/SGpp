/*
 * LejaPointDistribution.hpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_LEJAPOINTDISTRIBUTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_LEJAPOINTDISTRIBUTION_HPP_

#include "AbstractPointDistribution.hpp"
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/SingleFunction.hpp>
#include <vector>
#include <functional>
#include <cmath>

namespace SGPP {
namespace combigrid {

// for the test cases
void calc_leja_points(std::vector<float_t>& sortedPoints,
		std::vector<float_t> &points, int number, float_t lower_bound,
		float_t upper_bound, std::function<float_t(float_t)> weight_func);

/**
 * Provides Leja points (which are nested, i. e. the set of n leja points is a subset of the set of n+1 leja points)
 */
class LejaPointDistribution: public AbstractPointDistribution {
	std::vector<SGPP::float_t> points;
	std::vector<SGPP::float_t> sortedPoints;
	SingleFunction weightFunction;

	SGPP::float_t startingPoint;

	SGPP::float_t calcStartingPoint();
	// normal distribution to weight the weight function
	std::function<SGPP::float_t (SGPP::float_t)> normalDistribution;

public:
	// TODO: add default constructor with constant weight function and precomputed values?
	LejaPointDistribution(SingleFunction weightFunction = SingleFunction(constantFunction<float_t>(static_cast<float_t>(1.0))));
	virtual ~LejaPointDistribution();

	virtual SGPP::float_t compute(size_t numPoints, size_t j);
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_LEJAPOINTDISTRIBUTION_HPP_ */

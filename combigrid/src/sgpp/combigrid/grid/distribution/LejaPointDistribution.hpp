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

namespace sgpp{
namespace combigrid {

// for the test cases
void calc_leja_points(std::vector<double>& sortedPoints,
		std::vector<double> &points, int number, double lower_bound,
		double upper_bound, std::function<double(double)> weight_func);

/**
 * Provides Leja points (which are nested, i. e. the set of n leja points is a subset of the set of n+1 leja points)
 */
class LejaPointDistribution: public AbstractPointDistribution {
	std::vector<double> points;
	std::vector<double> sortedPoints;
	SingleFunction weightFunction;

	double startingPoint;

	double calcStartingPoint();
	// normal distribution to weight the weight function
	std::function<double (double)> normalDistribution;

public:
	// TODO: add default constructor with constant weight function and precomputed values?
	LejaPointDistribution(SingleFunction weightFunction = SingleFunction(constantFunction<double>(static_cast<double>(1.0))));
	virtual ~LejaPointDistribution();

	virtual double compute(size_t numPoints, size_t j);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_LEJAPOINTDISTRIBUTION_HPP_ */

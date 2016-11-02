// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_LEJAPOINTDISTRIBUTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_LEJAPOINTDISTRIBUTION_HPP_

#include <sgpp/combigrid/SingleFunction.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/distribution/AbstractPointDistribution.hpp>

#include <cmath>
#include <functional>
#include <vector>

namespace sgpp {
namespace combigrid {

// used only internally and for the test cases
void calc_leja_points(std::vector<double>& sortedPoints, std::vector<double>& points, int number,
                      double lower_bound, double upper_bound,
                      std::function<double(double)> weight_func);

/**
 * Provides Leja points (which are nested, i. e. the set of n leja points is a subset of the set of
 * n+1 leja points). Also accepts a weight function to compute weighted Leja points.
 */
class LejaPointDistribution : public AbstractPointDistribution {
  std::vector<double> points;
  std::vector<double> sortedPoints;
  SingleFunction weightFunction;

  double startingPoint;

  double calcStartingPoint();
  // normal distribution to weight the weight function
  std::function<double(double)> normalDistribution;

 public:
  // TODO(holzmudd): add default constructor with constant weight function and precomputed values?
  LejaPointDistribution(SingleFunction weightFunction =
                            SingleFunction(constantFunction<double>(static_cast<double>(1.0))));
  virtual ~LejaPointDistribution();

  virtual double compute(size_t numPoints, size_t j);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_DISTRIBUTION_LEJAPOINTDISTRIBUTION_HPP_ */

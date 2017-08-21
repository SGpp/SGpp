// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/optimization/MixedOptimizer.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * calculates Leja points
 * @param sortedPoints current Points
 * @param points current Points, new points will be added to this list
 * @param number number of next leja points, new points will be addes to this list
 * @param lower_bound lower bound for the range of the points
 * @param upper_bound upper bound for the range of the points
 * @param weight_func the weight function
 */
void LejaPointDistribution::calc_leja_points(std::vector<double>& sortedPoints,
                                             std::vector<double>& points, int number,
                                             double lower_bound, double upper_bound,
                                             std::function<double(double)> weight_func,
                                             double epsilon) {
  // calculates the next NUMBER leja points
  for (int i = 0; i < number; ++i) {
    // find the local optimum
    double x_opt = 0.0;
    double y_min = 0.0;

    // define remainder polynomial based on the current leja sequence
    // stored in sortedPoints
    std::function<double(double)> leja_func = [sortedPoints, weight_func](double x) {
      double prod = 1;
      for (size_t i = 0; i < sortedPoints.size(); ++i) {
        prod *= std::abs(x - sortedPoints.at(i));
      }
      return -prod * weight_func(x);
    };

    // optimize each interval of the remainder polynomial
    double x_lower = 0.0;
    double x_upper = lower_bound;
    for (size_t j = 0; j <= sortedPoints.size(); ++j) {
      // update the interval
      x_lower = x_upper;
      x_upper = (j < sortedPoints.size()) ? sortedPoints[j] : upper_bound;

      double x_val = (x_upper + x_lower) / 2.;
      double y_val = 0.0;

      // optimize the remainder polynomial if the current patch is wide enough
      if (std::abs(x_lower - x_upper) > 1e-10) {
        auto myLejaFunc = SingleFunction(leja_func);
        auto result = MixedOptimizer(myLejaFunc)
                          .minimize(OptimizationGuess::initial(x_lower, x_upper, myLejaFunc));
        x_val = result.b;
        y_val = result.fb;
      }

      // check if the maximum of the current patch is larger then
      // the maximum of the previous ones
      if (y_val < y_min) {
        x_opt = x_val;
        y_min = y_val;
      }
    }

    // add the current optimum to the sequence
    if (x_opt < 0 || x_opt > 1) {
      std::cout << "stop right here" << std::endl;
    }
    points.push_back(x_opt);

    // insertion sort
    for (size_t j = 0; j < sortedPoints.size(); ++j) {
      if (sortedPoints.at(j) > x_opt) {
        sortedPoints.insert(sortedPoints.begin() + j, x_opt);
        break;
      }
    }

    if (sortedPoints.at(sortedPoints.size() - 1) < x_opt) {
      sortedPoints.push_back(x_opt);
    }
  }
}

/**
 * Calculates the Starting Point by weighting the weight function with a wide normal distribution
 * and searching via optimizer for the maximum
 */
double LejaPointDistribution::calcStartingPoint(double epsilon) {
  // optimize it
  double x_val = 0.5;
  auto result =
      MixedOptimizer(weightFunction).minimize(OptimizationGuess::initial(0.0, 1.0, weightFunction));
  x_val = result.b;
  return x_val;
}

LejaPointDistribution::LejaPointDistribution()
    : weightFunction(SingleFunction(constantFunction<double>(static_cast<double>(1.0)))),
      startingPoint(0.5) {
  // TODO(holzmudd): add precomputed points?
  points.push_back(this->startingPoint);
  sortedPoints.push_back(this->startingPoint);
}

LejaPointDistribution::LejaPointDistribution(SingleFunction weightFunction)
    : weightFunction(weightFunction), startingPoint() {
  startingPoint = calcStartingPoint();
  points.push_back(this->startingPoint);
  sortedPoints.push_back(this->startingPoint);
}

LejaPointDistribution::~LejaPointDistribution() {}

double LejaPointDistribution::compute(size_t numPoints, size_t j) {
  if (points.size() <= j) {
    calc_leja_points(sortedPoints, points, static_cast<int>(j + 1 - points.size()), 0.0, 1.0,
                     weightFunction.getStdFunction());
  }

  return points[j];
}

} /* namespace combigrid */
} /* namespace sgpp*/

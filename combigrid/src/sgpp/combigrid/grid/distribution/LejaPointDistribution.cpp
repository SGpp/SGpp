// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>

#ifdef USE_DLIB
#include <dlib/optimization.h>
#endif

#include <functional>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

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
#ifdef USE_DLIB
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
        try {
          y_val = dlib::find_min_single_variable(leja_func, x_val, x_lower, x_upper, epsilon, 500);
        } catch (dlib::optimize_single_variable_failure& e) {
        }
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

#else
  std::cerr << "Error in sgpp::combigrid::calc_leja_points: "
            << "SG++ was compiled without Dlib support!\n";
#endif /* USE_DLIB */
}

/**
 * Calculates the Starting Point by weighting the weight function with a wide normal distribution
 * and searching via optimizer for the maximum
 */
double LejaPointDistribution::calcStartingPoint(double epsilon) {
#ifdef USE_DLIB
  // weight the weight function with the normal distribution
  std::function<double(double)> w = [this](double x) {
    const double factor = 0.2;
    double evalNormal = std::exp(-(factor * (x - 0.5)) * (factor * (x - 0.5)));
    return -(evalNormal * this->weightFunction(x));
  };

  // optimize it
  double x_val = 0.5;
  try {
    dlib::find_min_single_variable(w, x_val, 0.0, 1.0, epsilon, 100, 0.5);
  } catch (dlib::optimize_single_variable_failure& e) {
  }
  return x_val;
#else
  std::cerr << "Error in sgpp::combigrid::calcPointDistribution: "
            << "SG++ was compiled without Dlib support!\n";
  return -1.0;
#endif /* USE_DLIB */
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
                     weightFunction.getLambdaExpression());
  }

  return points[j];
}

} /* namespace combigrid */
} /* namespace sgpp*/

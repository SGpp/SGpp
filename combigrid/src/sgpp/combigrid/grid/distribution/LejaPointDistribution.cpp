// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/optimization/Optimization.hpp>
#include <functional>
#include <vector>
#include <algorithm>
#include <cmath>

namespace sgpp {
namespace combigrid {

const double epsilon = 0.00001;

class leja_f : public optimize::func_base {
 public:
  explicit leja_f(std::function<double(double)> f) { this->func = f; }
  virtual ~leja_f() {}

  double operator()(double x) { return func(x); }

 private:
  std::function<double(double)> func;
};

/**
 * calculates Leja points
 * @param sortedPoints current Points
 * @param points current Points, new points will be added to this list
 * @param number number of next leja points, new points will be addes to this list
 * @param lower_bound lower bound for the range of the points
 * @param upper_bound upper bound for the range of the points
 * @param weight_func the weight function
 */
void calc_leja_points(std::vector<double>& sortedPoints, std::vector<double>& points, int number,
                      double lower_bound, double upper_bound,
                      std::function<double(double)> weight_func) {
  // calculates the next NUMBER leja points
  for (int i = 0; i < number; ++i) {
    std::vector<double> tries;
    std::vector<double> values;

    double low, up;
    for (size_t j = 0; j <= sortedPoints.size(); ++j) {
      low = (j == 0) ? lower_bound : sortedPoints.at(j - 1);
      up = (j == sortedPoints.size()) ? upper_bound : sortedPoints.at(j);

      std::function<double(double)> leja_func = [sortedPoints, weight_func](double x) {
        double prod = 1;
        for (size_t i = 0; i < sortedPoints.size(); ++i) {
          prod *= std::abs(x - sortedPoints.at(i));
        }
        return (-1) * prod * weight_func(x);
      };

      double x_val = (low + up) / 2;

      leja_f f(leja_func);
      double y_val = optimize::local_min(low, up, epsilon, f, x_val);

      tries.push_back(x_val);
      values.push_back(y_val);
    }

    double min_val = 100;
    size_t idx = 0;
    for (size_t j = 0; j < values.size(); ++j) {
      if (values.at(j) < min_val) {
        min_val = values.at(j);
        idx = j;
      }
    }

    double next_point = tries.at(idx);
    points.push_back(next_point);
    // insertion sort
    for (size_t j = 0; j < sortedPoints.size(); ++j) {
      if (sortedPoints.at(j) > next_point) {
        sortedPoints.insert(sortedPoints.begin() + j, next_point);
        break;
      }
    }
    if (sortedPoints.at(sortedPoints.size() - 1) < next_point) {
      sortedPoints.push_back(next_point);
    }
  }
}

/**
 * Calculates the Starting Point by weighting the weight function with a wide normal distribution
 * and searching via optimizer for the maximum
 */
double LejaPointDistribution::calcStartingPoint() {
  // weight the weight function with the normal distribution
  std::function<double(double)> w = [this](double x) {
    const double factor = 0.2;
    double evalNormal = std::exp(-(factor * (x - 0.5)) * (factor * (x - 0.5)));
    return -(evalNormal * this->weightFunction(x));
  };

  // optimize it
  leja_f f(w);
  // starting x value
  double x_val = 0.5;
  optimize::local_min(0.0, 1.0, epsilon, f, x_val);

  return x_val;
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

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/Distribution.hpp>

#include <iostream>
#include <random>

namespace sgpp {
namespace base {

/**
 */
class DistributionNormal : public Distribution {
 public:
  /**
   * Constructor
   */
  DistributionNormal(double mean, double stddev)
      : Distribution(), mean(mean), stddev(stddev), dist(mean, stddev) {}

  /**
   * Destructor
   */
  virtual ~DistributionNormal() {}

  /**
   *
   */
  double sample() { return dist(gen); }

  /**
   *
   */
  double eval(double x) {
    return 1.0 / (sqrt(2 * M_PI) * stddev) * exp(-(x - mean) * (x - mean) / (2 * stddev * stddev));
  }

  /**
   * Heuristical bounds.
   * For these bounds and quadrature order 20, for the objective function f(x)=1 mean and variance
   * are calculated with an error of machine precision => the normal distirbution is truncated s.t.
   * the missing part is not numerically relevant
   */
  sgpp::base::DataVector getBounds() {
    sgpp::base::DataVector bounds(2);
    bounds[0] = mean - 20 * stddev * stddev;
    bounds[1] = mean + 20 * stddev * stddev;
    return bounds;
  }

 private:
  double mean;
  double stddev;
  std::normal_distribution<double> dist;
};

}  // namespace base
}  // namespace sgpp

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
  DistributionNormal(double mean, double stddev, double l, double h)
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
    return 1 / (sqrt(2 * M_PI) * stddev) * exp(-(x - mean) * (x - mean) / (2 * stddev * stddev));
  }

  /**
   * Heuristical bounds.
   * For really accurate calculatins they must be expanded until less than 10^(-16) of the
   * distributions volume is outside the bounds
   */
  void getBounds(double& l, double& r) {
    l = mean - 10 * stddev * stddev;
    r = mean + 10 * stddev * stddev;
  }

 private:
  double mean;
  double stddev;
  std::normal_distribution<double> dist;
};

}  // namespace base
}  // namespace sgpp

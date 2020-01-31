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
 * Truncated normal distribution.
 * Only accepts samples within [lower,upper]
 */
class DistributionTruncNormal : public Distribution {
 public:
  /**
   * Constructor
   */
  DistributionTruncNormal(double mean, double stddev, double lower, double upper)
      : Distribution(),
        mean(mean),
        stddev(stddev),
        lower(lower),
        upper(upper),
        dist(mean, stddev) {}

  /**
   * Destructor
   */
  virtual ~DistributionTruncNormal() {}

  /**
   *
   */
  double sample() {
    while (true) {
      double number = dist(gen);
      if (number >= lower && number <= upper) return number;
    }
  }

  /**
   *
   */
  double eval(double x) {
    if (x < lower) {
      std::cout << "DistributionTruncNormal: argument not in specified interval!\n";
      return lower;
    } else if (x > upper) {
      std::cout << "DistributionTruncNormal: argument not in specified interval!\n";
      return upper;
    } else {
      return 1.0 / (sqrt(2 * M_PI) * stddev) *
             exp(-(x - mean) * (x - mean) / (2 * stddev * stddev));
    }
  }

  /**
   *
   */
  sgpp::base::DataVector getBounds() {
    sgpp::base::DataVector bounds(2);
    bounds[0] = lower;
    bounds[1] = upper;
    return bounds;
  }

  sgpp::base::DistributionType getType() { return sgpp::base::DistributionType::TruncNormal; }

  sgpp::base::DataVector getCharacteristics() {
    sgpp::base::DataVector characteristics(2);
    characteristics[0] = mean;
    characteristics[1] = stddev;
    return characteristics;
  }

 private:
  // mean, often called mu
  double mean;
  // standard deviation, often called sigma
  // Note: This is sigma, not sigma^2!
  double stddev;
  // lower bound
  double lower;
  // upper bound
  double upper;
  std::normal_distribution<double> dist;
};  // namespace base

}  // namespace base
}  // namespace sgpp

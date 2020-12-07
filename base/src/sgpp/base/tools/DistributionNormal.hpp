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
      : Distribution(), mean(mean), stddev(stddev), dist(mean, stddev) {
        // std::cout << "DistributionNormal: Overwritten Bounds!\n";
      }

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
   * According to Wikipedia (https://de.wikipedia.org/wiki/Normalverteilung#Streuintervalle)
   * inside the interval [mean-z*sigma,mean+z*sigma] lie 2*Phi(z)-1 of the mass of the normal
   * distribution.
   * For z = 7 less than 10^(-9) are outside [mean - 7*sigma, mean + 7*sigma]
   * For z = 8 less than 10^(-14) are outside [mean - 8*sigma, mean + 8*sigma]
   * For z = 9 less than 10^(-18) are outside [mean - 9*sigma, mean + 9*sigma]
   */
  sgpp::base::DataVector getBounds() {
    sgpp::base::DataVector bounds(2);
    bounds[0] = mean - 9 * stddev;
    bounds[1] = mean + 9 * stddev;

    // used these for the borehole example
    // bounds[0] = 0.05;
    // bounds[1] = 0.15;
    return bounds;
  }

  sgpp::base::DistributionType getType() { return sgpp::base::DistributionType::Normal; }

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
  std::normal_distribution<double> dist;
};

}  // namespace base
}  // namespace sgpp

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <sgpp/base/tools/Distribution.hpp>
#include <string>

namespace sgpp {
namespace base {

/**
 * Disclaimer: This implementation models the beta distribution used for Polynomial Chaos Expansion
 * and should not be confused with the textbook definition of a beta distribution
 */
class DistributionBeta : public Distribution {
 public:
  /**
   * Constructor
   */
  explicit DistributionBeta(double alpha, double beta)
      : Distribution(), alpha(alpha), beta(beta), dist1(alpha + 1), dist2(beta + 1) {}

  /**
   * Destructor
   */
  virtual ~DistributionBeta() {}

  /**
   *
   */
  double sample() {
    auto num1 = dist1(gen);
    auto num2 = dist2(gen);
    return (num1 / (num1 + num2)) * 2 - 1;
  }

  double eval(double x) {
    return (std::pow(1 - x, alpha) * std::pow(1 + x, beta)) /
           (std::pow(2, alpha + beta + 1) *
            ((std::tgamma(alpha + 1) * std::tgamma(beta + 1)) / std::tgamma(alpha + beta + 2)));
  }

  sgpp::base::DataVector getBounds() {
    sgpp::base::DataVector bounds{-1, 1};
    return bounds;
  }

  sgpp::base::DistributionType getType() { return sgpp::base::DistributionType::Beta; }

  sgpp::base::DataVector getCharacteristics() { return sgpp::base::DataVector{alpha, beta}; }

 private:
  double alpha;
  double beta;
  std::gamma_distribution<double> dist1;
  std::gamma_distribution<double> dist2;
};
}  // namespace base
}  // namespace sgpp

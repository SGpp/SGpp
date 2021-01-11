// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <cmath>
#include <iostream>
#include <random>
#include <sgpp/base/tools/Distribution.hpp>
#include <string>

namespace sgpp {
namespace base {

/**
 */
class DistributionTruncGamma : public Distribution {
 public:
  /**
   * Constructor
   */
  explicit DistributionTruncGamma(double alpha = 0, double r = 20)
      : Distribution(), alpha(alpha), r(r), dist(alpha) {}

  /**
   * Destructor
   */
  virtual ~DistributionTruncGamma() {}

  /**
   *
   */
  double sample() {
    while (true) {
      auto num = dist(gen);
      if (num <= r) {
        return num;
      }
    }
  }
  double eval(double x) { return (std::pow(x, alpha) * std::exp(-x)) / std::tgamma(alpha + 1); }

  sgpp::base::DataVector getBounds() {
    sgpp::base::DataVector bounds{0, r};
    return bounds;
  }

  sgpp::base::DistributionType getType() { return sgpp::base::DistributionType::TruncGamma; }

  sgpp::base::DataVector getCharacteristics() { return sgpp::base::DataVector{alpha, r}; }

 private:
  double alpha;
  double r;
  std::gamma_distribution<double> dist;
};
}  // namespace base
}  // namespace sgpp

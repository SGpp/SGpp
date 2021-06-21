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
class DistributionTruncExponential : public Distribution {
 public:
  /**
   * Constructor
   */
  explicit DistributionTruncExponential(double r = 20, double lambda = 1)
      : Distribution(), lambda(lambda), r(r), dist(lambda) {}

  /**
   * Destructor
   */
  virtual ~DistributionTruncExponential() {}

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
  double eval(double x) { return lambda * std::exp(-x * lambda); }

  sgpp::base::DataVector getBounds() {
    sgpp::base::DataVector bounds{0, r};
    return bounds;
  }

  sgpp::base::DistributionType getType() { return sgpp::base::DistributionType::TruncExponential; }

  sgpp::base::DataVector getCharacteristics() { return sgpp::base::DataVector{r, lambda}; }

 private:
  double lambda;
  double r;
  std::exponential_distribution<double> dist;
};
}  // namespace base
}  // namespace sgpp

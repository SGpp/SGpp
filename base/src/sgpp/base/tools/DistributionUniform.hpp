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
class DistributionUniform : public Distribution {
 public:
  /**
   * Constructor
   */
  explicit DistributionUniform(double l = 0, double r = 1)
      : Distribution(), l(l), r(r), dist(l, r) {}

  /**
   * Destructor
   */
  virtual ~DistributionUniform() {}

  /**
   *
   */
  double sample() { return dist(gen); }

  double eval(double x) { return 1.0 / (r - l); }

  void getBounds(double& l, double& r) {
    l = this->l;
    r = this->r;
  }

 private:
  double l;
  double r;
  std::uniform_real_distribution<double> dist;
};
}  // namespace base
}  // namespace sgpp

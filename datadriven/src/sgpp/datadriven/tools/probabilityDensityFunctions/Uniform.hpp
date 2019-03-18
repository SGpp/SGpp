// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/tools/probabilityDensityFunctions/ProbabilityDensityFunction.hpp>

#include <iostream>
#include <random>

namespace sgpp {
namespace datadriven {

/**
 */
class Uniform : public ProbabilityDensityFunction {
 public:
  /**
   * Constructor
   */
  explicit Uniform(double l = 0, double h = 1)
      : ProbabilityDensityFunction(), l(l), h(h), dist(l, h) {}

  /**
   * Destructor
   */
  virtual ~Uniform() {}

  /**
   *
   */
  double sample() { return dist(gen); }

 private:
  double l;
  double h;
  std::uniform_real_distribution<double> dist;
};
}  // namespace datadriven
}  // namespace sgpp

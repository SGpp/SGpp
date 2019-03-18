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
class Normal : public ProbabilityDensityFunction {
 public:
  /**
   * Constructor
   */
  Normal(double mean, double stddev)
      : ProbabilityDensityFunction(), mean(mean), stddev(stddev), dist(mean, stddev) {}

  /**
   * Destructor
   */
  virtual ~Normal() {}

  /**
   *
   */
  double sample() { return dist(gen); }

 private:
  double mean;
  double stddev;
  std::normal_distribution<double> dist;
};

}  // namespace datadriven
}  // namespace sgpp

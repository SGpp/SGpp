// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>

#include <chrono>
#include <iostream>
#include <random>

namespace sgpp {
namespace datadriven {

/**
 * stores a sparse grid not a knot B-spline interpolant in the framework of a respsonse surface
 */
class ProbabilityDensityFunction {
 public:
  /**
   * Constructor
   */
  ProbabilityDensityFunction() {
    // set seed
    gen.seed(std::chrono::system_clock::now().time_since_epoch().count());
  }

  /**
   * Destructor
   */
  virtual ~ProbabilityDensityFunction() {}

  /**
   *
   */
  virtual double sample() = 0;

  sgpp::base::DataVector samples(size_t num);

 protected:
  std::default_random_engine gen;
};

}  // namespace datadriven
}  // namespace sgpp

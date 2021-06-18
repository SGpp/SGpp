// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <chrono>
#include <iostream>
#include <random>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <string>

namespace sgpp {
namespace base {

/**
 * enum to address different pdfs in a standardized way
 *
 */
enum class DistributionType {
  Uniform,      //  0
  Normal,       //  1
  Lognormal,    //  2
  TruncNormal,  // 3
  Beta,
  TruncExponential,
  TruncGamma
};

/**
 * stores a sparse grid not a knot B-spline interpolant in the framework of a respsonse surface
 */
class Distribution {
 public:
  /**
   * Constructor
   *
   * @param seed  if a seed should be set, ii.e. for precalculating and reusing grids set one here,
   * otherwise a pseudo random seed is set automatically
   */
  Distribution(typename std::chrono::system_clock::duration::rep seed = 777) {
    // set seed
    if (seed == 777) {
      seed = std::chrono::system_clock::now().time_since_epoch().count();
    }
    gen.seed(static_cast<decltype(gen)::result_type>(seed));
  }

  /**
   * Destructor
   */
  virtual ~Distribution() {}

  /**
   *
   */
  virtual double sample() = 0;

  /**
   *
   */
  virtual double eval(double x) = 0;

  sgpp::base::DataVector samples(size_t num);

  virtual sgpp::base::DataVector getBounds() = 0;

  virtual sgpp::base::DistributionType getType() = 0;

  /**
   * return all relevant characteristic values (e.g. mean and standarddeviation for normal
   * distribution)
   */
  virtual sgpp::base::DataVector getCharacteristics() = 0;

 protected:
  std::default_random_engine gen;
};

}  // namespace base
}  // namespace sgpp

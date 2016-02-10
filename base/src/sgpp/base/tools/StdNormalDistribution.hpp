// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STDNORMALDISTRIBUTION_HPP
#define STDNORMALDISTRIBUTION_HPP

#include <sgpp/globaldef.hpp>

#include <cmath>


namespace SGPP {
namespace base {

/**
 * This provides an implementation of the standradized normal
 * distribution.
 *
 * Also a density function for the normal distribution is provided!
 */
class StdNormalDistribution {
 public:
  /**
   * Std-Constructor
   */
  StdNormalDistribution();

  /**
   * Std-Destructor
   */
  ~StdNormalDistribution();

  /**
   * Calculates the Cumulative Density values of the standard normal distribution
   * (expected values = 0.0, standard deviation = 1.0)
   *
   * taken from http://www.richelbilderbeek.nl/CppGetCumulativeDensityNormal.htm
   *
   * @param x the value for which the cumulative density value should be calculated
   */
  float_t getCumulativeDensity(const float_t x);

  /*
   * Calculates the Density values of the standard normal distribution
   * (expected values = 0.0, standard deviation = 1.0)
   *
   * taken from http://www.richelbilderbeek.nl/CppGetDensityNormal.htm
   *
   * @param x the value for which the density value should be calculated
   */
  float_t getDensity(const float_t x);

  /*
   * Calculates the Density values of the normal distribution
   *
   * taken from http://www.richelbilderbeek.nl/CppGetDensityNormal.htm
   *
   * @param x the value for which the density value should be calculated
   * @param mu the expected value of the normal distribution
   * @param sigma the standard deviation of the normal distribution
   */
  float_t getDensity(const float_t x, const float_t mu, const float_t sigma);


  /*
   * Calculates the Density values of the normal distribution (normed variant)
   *
   * @param x the value for which the density value should be calculated
   * @param mu the expected value of the normal distribution
   * @param sigma the standard deviation of the normal distribution
   */
  float_t getNormedDensity(const float_t x, const float_t mu,
                           const float_t sigma);
};

}  // namespace base
}  // namespace SGPP

#endif /* STDNORMALDISTRIBUTION_HPP */

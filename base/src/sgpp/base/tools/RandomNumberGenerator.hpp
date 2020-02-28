// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>

#include <random>

namespace sgpp {
namespace base {

/**
 * Singleton class for generating pseudo-random numbers
 * (wrapper around std::mt19937 from &lt;random&gt;).
 */
class RandomNumberGenerator {
 public:
  /// type of the seed
  typedef std::mt19937::result_type SeedType;

  /**
   * @return singleton instance
   */
  static RandomNumberGenerator& getInstance();

  /**
   * Generate a uniform pseudo-random number.
   *
   * @param a lower bound
   * @param b upper bound
   * @return  uniform pseudo-random number in \f$[a, b]\f$
   */
  double getUniformRN(double a = 0.0, double b = 1.0);

  /**
   * Fills vector with uniform pseudo-random numbers.
   *
   * @param[out]  vector  vector to be filled, has to have desired size beforehand
   * @param       a       lower bound
   * @param       b       upper bound
   */
  void getUniformRV(DataVector& vector, double a = 0.0, double b = 1.0);

  /**
   * Generate a uniform pseudo-random array index.
   *
   * @param size  size of the array
   * @return      discrete uniform pseudo-random number in
   *              \f$\{0, \dotsc, \text{\texttt{size}} - 1\}\f$
   */
  size_t getUniformIndexRN(size_t size);

  /**
   * Generate a Gaussian pseudo-random number.
   *
   * @param mean      mean of the Gaussian distribution
   * @param stdDev    standard deviation of the Gaussian distribution
   * @return          Gaussian pseudo-random number
   */
  double getGaussianRN(double mean = 0.0, double stdDev = 1.0);

  /**
   * Fills vector with Gaussian pseudo-random numbers.
   *
   * @param[out]  vector  vector to be filled, has to have desired size beforehand
   * @param       mean    mean of the Gaussian distribution
   * @param       stdDev  standard deviation of the Gaussian distribution
   */
  void getGaussianRV(DataVector& vector, double mean = 0.0, double stdDev = 1.0);

  /**
   * @return      seed
   */
  SeedType getSeed() const;

  /**
   * Reseeds with time-dependent seed.
   */
  void setSeed();

  /**
   * Reseeds.
   *
   * @param seed  seed to be used
   */
  void setSeed(SeedType seed);

 protected:
  /// random number generator
  std::mt19937 generator;
  /// seed
  SeedType seed;

 private:
  /**
   * Constructor, initializes with time-dependent seed.
   */
  RandomNumberGenerator();

  /**
   * Deleted copy constructor.
   */
  RandomNumberGenerator(const RandomNumberGenerator&) = delete;

  /**
   * Deleted assignment operator.
   */
  void operator=(const RandomNumberGenerator&) = delete;
};
}  // namespace base
}  // namespace sgpp

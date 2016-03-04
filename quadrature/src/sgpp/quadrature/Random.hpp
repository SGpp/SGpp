// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <sgpp/globaldef.hpp>

#include <random>

namespace sgpp {
namespace quadrature {

/**
 * Wraps the random generator to use. Ensures that it is is seeded correctly.
 */

class Random {
 public:
  /**
   * Manually seed the generator with a given seed.
   * @param seed_value the seed to use.
   */
  static void seed(std::uint64_t seed_value = std::mt19937_64::default_seed);

  /**
   * returns a random integer value in [0, RAND_MAX)
   */
  static std::uint64_t random_uint64();

  /**
   * returns a random double value in [0, 1)
   */
  static double random_double();

 protected:
  static bool is_seeded;

  static std::mt19937_64 gen;
  static std::uniform_int_distribution<std::uint64_t> distInt;
  static std::uniform_real_distribution<double> distReal;
};
}  // namespace quadrature
}  // namespace sgpp

#endif /* RANDOM_HPP */

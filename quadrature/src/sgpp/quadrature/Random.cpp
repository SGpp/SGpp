// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/quadrature/Random.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace quadrature {

bool Random::is_seeded = false;
std::mt19937_64 Random::gen = std::mt19937_64();
std::uniform_int_distribution<std::uint64_t> Random::distInt =
    std::uniform_int_distribution<std::uint64_t>(0, RAND_MAX);
std::uniform_real_distribution<double> Random::distReal =
    std::uniform_real_distribution<double>(0, 1);

void Random::seed(std::uint64_t seed_value) {
  gen.seed(seed_value);
  is_seeded = true;
}

std::uint64_t Random::random_uint64() {
  if (!is_seeded) {
    Random::seed(std::mt19937_64::default_seed);
  }

  return distInt(gen);
}

double Random::random_double() {
  if (!is_seeded) {
    Random::seed(std::mt19937_64::default_seed);
  }

  return distReal(gen);
}

}  // namespace quadrature
}  // namespace sgpp

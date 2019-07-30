// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/globaldef.hpp>

#include <ctime>

namespace sgpp {
namespace base {

RandomNumberGenerator::RandomNumberGenerator() { setSeed(); }

RandomNumberGenerator& RandomNumberGenerator::getInstance() {
  static RandomNumberGenerator rng;
  return rng;
}

double RandomNumberGenerator::getUniformRN(double a, double b) {
  std::uniform_real_distribution<double> distr(a, b);
  return distr(generator);
}

void RandomNumberGenerator::getUniformRV(base::DataVector& vector, double a, double b) {
  for (size_t i = 0; i < vector.getSize(); i++) {
    vector[i] = getUniformRN(a, b);
  }
}

size_t RandomNumberGenerator::getUniformIndexRN(size_t size) {
  std::uniform_int_distribution<size_t> distr(0, size - 1);
  return distr(generator);
}

double RandomNumberGenerator::getGaussianRN(double mean, double stdDev) {
  std::normal_distribution<double> distr(mean, stdDev);
  return distr(generator);
}

void RandomNumberGenerator::getGaussianRV(base::DataVector& vector, double mean, double stdDev) {
  for (size_t i = 0; i < vector.getSize(); i++) {
    vector[i] = getGaussianRN(mean, stdDev);
  }
}

RandomNumberGenerator::SeedType RandomNumberGenerator::getSeed() const { return seed; }

void RandomNumberGenerator::setSeed() { setSeed(static_cast<SeedType>(std::time(0))); }

void RandomNumberGenerator::setSeed(RandomNumberGenerator::SeedType seed) {
  this->seed = seed;
  generator.seed(seed);
}
}  // namespace base
}  // namespace sgpp

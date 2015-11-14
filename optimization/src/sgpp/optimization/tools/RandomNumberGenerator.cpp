// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <ctime>

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

namespace SGPP {
  namespace optimization {

    RandomNumberGenerator::RandomNumberGenerator() {
      setSeed();
    }

    float_t RandomNumberGenerator::getUniformRN(float_t a, float_t b) {
      std::uniform_real_distribution<float_t> distr(a, b);
      return distr(generator);
    }

    size_t RandomNumberGenerator::getUniformIndexRN(size_t size) {
      std::uniform_int_distribution<size_t> distr(0, size - 1);
      return distr(generator);
    }

    float_t RandomNumberGenerator::getGaussianRN(float_t stdDev,
        float_t mean) {
      std::normal_distribution<float_t> distr(mean, stdDev);
      return distr(generator);
    }

    RandomNumberGenerator::SeedType RandomNumberGenerator::getSeed() const {
      return seed;
    }

    void RandomNumberGenerator::setSeed() {
      setSeed(static_cast<SeedType>(std::time(0)));
    }

    void RandomNumberGenerator::setSeed(RandomNumberGenerator::SeedType seed) {
      this->seed = seed;
      generator.seed(seed);
    }

  }
}

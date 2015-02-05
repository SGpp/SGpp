// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <ctime>

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/tools/RNG.hpp>

namespace SGPP {
  namespace optimization {
    namespace tools {

      RNG rng;

      RNG::RNG() {
        setSeed();
      }

      float_t RNG::getUniformRN(float_t a, float_t b) {
        std::uniform_real_distribution<float_t> distr(a, b);
        return distr(generator);
      }

      size_t RNG::getUniformIndexRN(size_t size) {
        std::uniform_int_distribution<size_t> distr(0, size - 1);
        return distr(generator);
      }

      float_t RNG::getGaussianRN(float_t stdDev, float_t mean) {
        std::normal_distribution<float_t> distr(mean, stdDev);
        return distr(generator);
      }

      RNG::SeedType RNG::getSeed() {
        return seed;
      }

      void RNG::setSeed() {
        setSeed(static_cast<SeedType>(std::time(0)));
      }

      void RNG::setSeed(RNG::SeedType seed) {
        this->seed = seed;
        generator.seed(seed);
      }

    }
  }
}

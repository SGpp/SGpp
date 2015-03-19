// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/StdNormalDistribution.hpp>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    StdNormalDistribution::StdNormalDistribution() {
    }

    StdNormalDistribution::~StdNormalDistribution() {
    }

    float_t StdNormalDistribution::getCumulativeDensity(const float_t x) {
      const float_t c0 = 0.2316419;
      const float_t c1 = 1.330274429;
      const float_t c2 = 1.821255978;
      const float_t c3 = 1.781477937;
      const float_t c4 = 0.356563782;
      const float_t c5 = 0.319381530;
      const float_t c6 = 0.398942280401;

      const float_t negative = (x < 0 ? 1.0 : 0.0);
      const float_t xPos = (x < 0.0 ? -x : x);
      const float_t k = 1.0 / ( 1.0 + (c0 * xPos));
      const float_t y1 = (((((((c1 * k - c2) * k) + c3) * k) - c4) * k) + c5) * k;
      const float_t y2 = 1.0 - (c6 * std::exp(-0.5 * xPos * xPos) * y1);

      return ((1.0 - negative) * y2) + (negative * (1.0 - y2));
    }


    float_t StdNormalDistribution::getDensity(const float_t x) {
      const float_t mean = 0.0;
      const float_t stddev = 1.0;

      const float_t firstTerm = 1.0 / (stddev * std::sqrt(2.0 * M_PI));
      const float_t secondTerm = -( (x - mean) * (x - mean) / (2.0 * stddev * stddev) );
      const float_t result = firstTerm * std::exp(secondTerm);

      return result;
    }

    float_t StdNormalDistribution::getDensity(const float_t x, const float_t mu, const float_t sigma) {
      const float_t mean = mu;
      const float_t stddev = sigma;

      const float_t firstTerm = 1.0 / (stddev * std::sqrt(2.0 * M_PI));
      const float_t secondTerm = -( (x - mean) * (x - mean) / (2.0 * stddev * stddev) );
      const float_t result = firstTerm * std::exp(secondTerm);

      return result;
    }

    float_t StdNormalDistribution::getNormedDensity(const float_t x, const float_t mu, const float_t sigma) {
      const float_t mean = mu;
      const float_t stddev = sigma;

      const float_t secondTerm = -( (x - mean) * (x - mean) / (2.0 * stddev * stddev) );
      const float_t result = std::exp(secondTerm);

      return result;
    }

  }
}
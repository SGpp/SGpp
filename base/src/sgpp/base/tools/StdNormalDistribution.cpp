/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

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

    double StdNormalDistribution::getCumulativeDensity(const double x) {
      const double c0 = 0.2316419;
      const double c1 = 1.330274429;
      const double c2 = 1.821255978;
      const double c3 = 1.781477937;
      const double c4 = 0.356563782;
      const double c5 = 0.319381530;
      const double c6 = 0.398942280401;

      const double negative = (x < 0 ? 1.0 : 0.0);
      const double xPos = (x < 0.0 ? -x : x);
      const double k = 1.0 / ( 1.0 + (c0 * xPos));
      const double y1 = (((((((c1 * k - c2) * k) + c3) * k) - c4) * k) + c5) * k;
      const double y2 = 1.0 - (c6 * std::exp(-0.5 * xPos * xPos) * y1);

      return ((1.0 - negative) * y2) + (negative * (1.0 - y2));
    }


    double StdNormalDistribution::getDensity(const double x) {
      const double mean = 0.0;
      const double stddev = 1.0;

      const double firstTerm = 1.0 / (stddev * std::sqrt(2.0 * M_PI));
      const double secondTerm = -( (x - mean) * (x - mean) / (2.0 * stddev * stddev) );
      const double result = firstTerm * std::exp(secondTerm);

      return result;
    }

    double StdNormalDistribution::getDensity(const double x, const double mu, const double sigma) {
      const double mean = mu;
      const double stddev = sigma;

      const double firstTerm = 1.0 / (stddev * std::sqrt(2.0 * M_PI));
      const double secondTerm = -( (x - mean) * (x - mean) / (2.0 * stddev * stddev) );
      const double result = firstTerm * std::exp(secondTerm);

      return result;
    }

    double StdNormalDistribution::getNormedDensity(const double x, const double mu, const double sigma) {
      const double mean = mu;
      const double stddev = sigma;

      const double secondTerm = -( (x - mean) * (x - mean) / (2.0 * stddev * stddev) );
      const double result = std::exp(secondTerm);

      return result;
    }

  }
}

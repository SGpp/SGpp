// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "HaltonSampleGenerator.hpp"
#include <sgpp/quadrature/Random.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>

using namespace SGPP::base;

namespace SGPP {
  namespace quadrature {

    HaltonSampleGenerator::HaltonSampleGenerator(size_t dimensions) : SampleGenerator(dimensions) {
      int basePrimes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};

      base_vektor = new int[dimensions];

      f_vektor = new float_t[dimensions];
      i_vektor = new float_t[dimensions];
      result_vektor = new float_t[dimensions];

      for (size_t i = 0; i < dimensions; i++) {
        base_vektor[i] = basePrimes[Random::random() % 15];
        f_vektor[i] = 1. / ((float_t)base_vektor[i]);

        index = 1;
        result_vektor[i] = 0.;
      }
    }

    void HaltonSampleGenerator::getSample(SGPP::base::DataVector& dv) {

      for (size_t i = 0; i < dimensions; i++) {

        result_vektor[i] = 0.;
        f_vektor[i] = 1. / ((float_t)base_vektor[i]);
        i_vektor[i] = index;

        while (i_vektor[i] > 0) {
          result_vektor[i] = result_vektor[i] + f_vektor[i] * ((float_t)((int)i_vektor[i] % base_vektor[i]));
          i_vektor[i] = floor(((float_t)i_vektor[i]) / ((float_t)base_vektor[i]));
          f_vektor[i] = f_vektor[i] / ((float_t)base_vektor[i]);
        }

        dv[i] = result_vektor[i];
      }

      index++;
    }

  }
}

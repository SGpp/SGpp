// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/quadrature/sample/NaiveSampleGenerator.hpp>

#include <sgpp/quadrature/Random.hpp>

using namespace SGPP::base;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace quadrature {

    void NaiveSampleGenerator::getSample(SGPP::base::DataVector& dv) {
      // generate random sample with dimensionality corresponding to the
      // size of the given datavector (in 0 to 1)
      for (size_t i = 0; i < dv.getSize(); i++) {
        dv[i] = Random::random_double();
      }
    }

  }
}
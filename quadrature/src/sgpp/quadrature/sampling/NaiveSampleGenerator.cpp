// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "NaiveSampleGenerator.hpp"

#include <sgpp/quadrature/Random.hpp>
#include <sgpp/globaldef.hpp>

using namespace SGPP::base;

namespace SGPP {
  namespace quadrature {

    NaiveSampleGenerator::NaiveSampleGenerator(size_t dimension,
        std::uint64_t seed) :
      SampleGenerator(dimension, seed),
      uniformRealDist(0, 1) {
    }

    NaiveSampleGenerator::~NaiveSampleGenerator() {
    }

    void NaiveSampleGenerator::getSample(DataVector& sample) {
      // generate random sample with dimensionality corresponding to the
      // size of the given datavector (in 0 to 1)
      for (size_t i = 0; i < sample.getSize(); i++) {
        sample[i] = uniformRealDist(rng);
      }
    }

  }
}

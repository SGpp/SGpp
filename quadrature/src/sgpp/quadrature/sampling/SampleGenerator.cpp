// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/quadrature/sampling/SampleGenerator.hpp>

#include <sgpp/quadrature/Random.hpp>
#include <sgpp/globaldef.hpp>

using namespace SGPP::base;

namespace SGPP {
  namespace quadrature {

    SampleGenerator::SampleGenerator(size_t dimensions, int seed) :
      dimensions(dimensions), seed(seed) {
      Random::seed(seed);
    }

    SampleGenerator::~SampleGenerator() {
    }

    void SampleGenerator::getSamples(DataMatrix& samples) {

      // Number of columns has to correspond to the number of dimensions
      if (samples.getNcols() != dimensions)
        return;

      // generate one sample for every row of the given DataMatrix
      DataVector dv(dimensions);

      for (size_t i = 0; i < samples.getNrows(); i++) {
        getSample(dv);
        samples.setRow(i, dv);
      }
    }

    size_t SampleGenerator::getDimensions() {
      return dimensions;
    }

    void SampleGenerator::setDimensions(size_t dimensions) {
      this->dimensions = dimensions;
    }

  }
}

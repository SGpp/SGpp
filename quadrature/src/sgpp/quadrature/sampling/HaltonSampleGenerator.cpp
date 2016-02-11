// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/quadrature/sampling/HaltonSampleGenerator.hpp>
#include <sgpp/globaldef.hpp>

#include <random>
#include <cmath>

namespace SGPP {
namespace quadrature {

HaltonSampleGenerator::HaltonSampleGenerator(size_t dimensions,
    std::uint64_t seed) :
  SampleGenerator(dimensions, seed),
  index(1),
  baseVector(dimensions),
  iVector(dimensions),
  fVector(dimensions),
  resultVector(dimensions),
  distInt(0, 15) {
  size_t basePrimes[] =
  { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47 };

  for (size_t i = 0; i < dimensions; i++) {
    baseVector[i] = basePrimes[distInt(rng)];
    fVector[i] = 1. / static_cast<float_t>(baseVector[i]);
    resultVector[i] = 0.;
  }
}

HaltonSampleGenerator::~HaltonSampleGenerator() {
}

void HaltonSampleGenerator::getSample(SGPP::base::DataVector& dv) {
  for (size_t i = 0; i < dimensions; i++) {
    resultVector[i] = 0.;
    fVector[i] = 1. / static_cast<float_t>(baseVector[i]);
    iVector[i] = static_cast<float_t>(index);

    while (iVector[i] > 0) {
      resultVector[i] = resultVector[i]
                        + fVector[i]
                        * ((float_t) ((size_t) iVector[i] % baseVector[i]));
      iVector[i] = floor(
                     ((float_t) iVector[i]) / ((float_t) baseVector[i]));
      fVector[i] = fVector[i] / ((float_t) baseVector[i]);
    }

    dv[i] = resultVector[i];
  }

  index++;
}

}  // namespace quadrature
}  // namespace SGPP

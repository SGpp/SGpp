// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/quadrature/sampling/StratifiedSampleGenerator.hpp>
#include <sgpp/quadrature/Random.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>
#include <vector>

namespace sgpp {
namespace quadrature {

StratifiedSampleGenerator::StratifiedSampleGenerator(std::vector<size_t>& strataPerDimension,
                                                     uint64_t seed)
    : SampleGenerator(strataPerDimension.size(), seed),
      numberOfStrata(strataPerDimension),
      currentStrata(strataPerDimension.size()),
      sizeOfStrata(strataPerDimension.size()),
      uniformRealDist(0, 1) {
  // set counter to the first strata for each dimension
  // compute size of strata per dimension
  for (size_t i = 0; i < dimensions; i++) {
    if (numberOfStrata[i] < 1) numberOfStrata[i] = 1;

    currentStrata[i] = 0;
    sizeOfStrata[i] = 1. / static_cast<double>(this->numberOfStrata[i]);
  }
}

StratifiedSampleGenerator::~StratifiedSampleGenerator() {}

void StratifiedSampleGenerator::getSample(sgpp::base::DataVector& dv) {
  // Check for correct dimension of the parameter vector
  if (dv.getSize() != dimensions) return;

  // Choose a random number inside the stratum selected for this dimension
  for (size_t i = 0; i < dimensions; i++) {
    dv[i] = (static_cast<double>(currentStrata[i]) + uniformRealDist(rng)) * sizeOfStrata[i];
  }

  // continue to the next stratum used for the next sample
  getNextStrata();
}

void StratifiedSampleGenerator::getNextStrata() {
  for (size_t i = 0; i < dimensions; i++) {
    // next stratum in this dimension available
    if (currentStrata[i] < numberOfStrata[i] - 1) {
      currentStrata[i]++;
      break;
    } else {
      // no more stratum in this dimension, start at the first stratum of this
      // dimension and go on with iterating over the next dimension
      currentStrata[i] = 0;
    }
  }
}

}  // namespace quadrature
}  // namespace sgpp

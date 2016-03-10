// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/quadrature/sampling/LatinHypercubeSampleGenerator.hpp>
#include <sgpp/quadrature/Random.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>
#include <list>
#include <algorithm>
#include <vector>

namespace sgpp {
namespace quadrature {

LatinHypercubeSampleGenerator::LatinHypercubeSampleGenerator(size_t dimensions,
                                                             size_t numberOfStrata,
                                                             std::uint64_t seed)
    : SampleGenerator(dimensions, seed),
      numberOfStrata(
          numberOfStrata),       // each dimension is divided in n strata to provide n sample points
      numberOfCurrentSample(1),  // index number of current sample [1, n]
      // equidistant split of [0,1] in n strata -> size of one stratum = 1 / n
      sizeOfStrata(1. / static_cast<double>(numberOfStrata)),
      uniformRealDist(0, 1) {
  for (size_t i = 0; i < dimensions; i++) {
    currentStrata.push_back(std::vector<size_t>());

    for (size_t j = 0; j < numberOfStrata; j++) {
      currentStrata[i].push_back(j);
    }
  }

  shuffleStrataSequence();
}

LatinHypercubeSampleGenerator::~LatinHypercubeSampleGenerator() {}

void LatinHypercubeSampleGenerator::getSample(sgpp::base::DataVector& sample) {
  // compute random value inside the current stratum selected from the shuffled strata sequence
  for (size_t i = 0; i < dimensions; i++) {
    sample[i] =
        (static_cast<double>(currentStrata[i][numberOfCurrentSample - 1]) + uniformRealDist(rng)) *
        sizeOfStrata;
  }

  // select next sample from strata sequence.
  // If one sequence is complete shuffle strata to get a new one.
  if (numberOfCurrentSample < numberOfStrata) {
    numberOfCurrentSample++;
  } else {
    numberOfCurrentSample = 0;
    shuffleStrataSequence();
  }
}

void LatinHypercubeSampleGenerator::shuffleStrataSequence() {
  for (size_t i = 0; i < dimensions; i++) {
    std::shuffle(currentStrata[i].begin(), currentStrata[i].end(), rng);
  }
}

}  // namespace quadrature
}  // namespace sgpp

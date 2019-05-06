// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LATINHYPERCUBESAMPLEGENERATOR_HPP
#define LATINHYPERCUBESAMPLEGENERATOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/quadrature/sampling/SampleGenerator.hpp>

#include <list>
#include <vector>
#include <random>

namespace sgpp {
namespace quadrature {

/**
 * The class NaiveSampleGenerator implements a simple MonteCarlo sample
 * generator. A sample is generated using the standard random number
 * generator from cmath and transforming the values to double range 0.0 to
 * 1.0.
 */
class LatinHypercubeSampleGenerator : public SampleGenerator {
 public:
  /**
   * Standard constructor
   *
   * @param dimensions number of dimensions used for sample generation
   * @param numberOfStrata number of strata
   * @param seed custom seed (defaults to default seed of mt19937_64)
   */

  LatinHypercubeSampleGenerator(size_t dimensions, size_t numberOfStrata,
                                std::uint64_t seed = std::mt19937_64::default_seed);

  /**
   * Destructor
   */
  virtual ~LatinHypercubeSampleGenerator();

  /**
   * This method generates one sample .
   * Implementation of the abstract Method getSample from SampelGenerator.
   *
   * @param sample DataVector storing the new generated sample vector.
   */

  void getSample(sgpp::base::DataVector& sample);

 private:
  /**
   * This method generates one sample .
   * Implementation of the abstract Method getSample from SampelGenerator.
   */

  void shuffleStrataSequence();

  //
  size_t numberOfStrata;

  //
  size_t numberOfCurrentSample;

  //
  double sizeOfStrata;

  //
  std::vector<std::vector<size_t> > currentStrata;

  //
  std::uniform_real_distribution<double> uniformRealDist;
};

}  // namespace quadrature
}  // namespace sgpp

#endif /* LATINHYPERCUBESAMPLEGENERATOR_HPP */

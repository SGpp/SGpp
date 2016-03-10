// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HALTONSAMPLEGENERATOR_HPP
#define HALTONSAMPLEGENERATOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/quadrature/sampling/SampleGenerator.hpp>

#include <vector>

namespace sgpp {
namespace quadrature {

/**
 *
 */
class HaltonSampleGenerator : public SampleGenerator {
 public:
  /**
   * Standard constructor
   *
   * @param dimension number of dimensions used for sample generation
   * @param seed custom seed (defaults to default seed of mt19937_64)
   */
  explicit HaltonSampleGenerator(size_t dimension,
                                 std::uint64_t seed = std::mt19937_64::default_seed);

  /**
   * Destructor
   */
  ~HaltonSampleGenerator();

  /**
   * This method generates one sample .
   * Implementation of the abstract Method getSample from SampelGenerator.
   *
   * @param sample DataVector storing the new generated sample vector.
   */
  virtual void getSample(sgpp::base::DataVector& sample);

 private:
  size_t index;
  std::vector<size_t> baseVector;
  std::vector<double> iVector;
  std::vector<double> fVector;
  std::vector<double> resultVector;
  //
  std::uniform_int_distribution<std::uint64_t> distInt;
};

}  // namespace quadrature
}  // namespace sgpp

#endif /* HALTONSAMPLEGENERATOR_HPP */

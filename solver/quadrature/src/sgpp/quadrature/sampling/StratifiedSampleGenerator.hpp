// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STRATIFIEDSAMPLEGENERATOR_HPP
#define STRATIFIEDSAMPLEGENERATOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/quadrature/sampling/SampleGenerator.hpp>

#include <random>
#include <vector>

namespace sgpp {
namespace quadrature {

/**
 * The class StratifiedSampleGenerator subdivides every dimension in a given
 * number of strata. For each strata one sample point is generated. In case
 * one sample has already been generated for every strata, the next requested
 * sample will be placed into the first strata.
 */
class StratifiedSampleGenerator : public SampleGenerator {
 public:
  /**
   * Standard constructor
   *
   * @param strataPerDimension array holding the number of strata used to
   * subdivide the specific dimension
   * @param seed custom seed (defaults to default seed of mt19937_64)
   */

  StratifiedSampleGenerator(std::vector<size_t>& strataPerDimension,
                            std::uint64_t seed = std::mt19937_64::default_seed);

  /**
   * Destructor
   */
  virtual ~StratifiedSampleGenerator();

  /**
   * This method generates one sample .
   * Implementation of the abstract Method getSample from SampelGenerator.
   *
   * @param sample DataVector storing the new generated sample vector.
   */

  void getSample(sgpp::base::DataVector& sample);

 private:
  // Array containing the number of strata per dimension
  std::vector<size_t> numberOfStrata;
  // Array containing the current strata number for every dimension
  std::vector<size_t> currentStrata;

  // Array containing the size of dimension i strata when dividing [0,1] into numberOfStrata[i]
  std::vector<double> sizeOfStrata;

  /**
   * This method computes in which strata the next sample should be generated.
   * Dimension after dimension each stratum is used to generate one sample point. As soon
   * as one dimension is completed the algorithm will start at the beginning of this dimension and
   * counts up the next dimension by 1.
   */
  void getNextStrata();

  //
  std::uniform_real_distribution<double> uniformRealDist;
};

}  // namespace quadrature
}  // namespace sgpp

#endif /* STRATIFIEDSAMPLEGENERATOR_HPP */

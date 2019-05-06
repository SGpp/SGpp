// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATUREMCADVANCED_HPP
#define OPERATIONQUADRATUREMCADVANCED_HPP

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/quadrature/sampling/SampleGenerator.hpp>

#include <vector>

namespace sgpp {
namespace quadrature {

/**
 * Typedef for general functions that can be passed to integration methods. Requires three
 * parameters. First, the dimensionality, then dim-many coordinates, and then further client data
 * for the function at hand.
 */
typedef double (*FUNC)(int, double*, void*);

/**
 * Quadrature on any sparse grid (that has OperationMultipleEval implemented)
 * using various Monte Carlo Methods (Advanced).
 */

class OperationQuadratureMCAdvanced : public sgpp::base::OperationQuadrature {
 public:
  /**
   * @brief Constructor of OperationQuadratureMCAdvanced, specifying a grid
   * object and the number of samples to use.
   *
   * @param grid Reference to the grid object
   * @param numberOfSamples Number of Monte Carlo samples
   * @param seed Custom seed (defaults to default seed of mt19937_64)
   */
  OperationQuadratureMCAdvanced(sgpp::base::Grid& grid, size_t numberOfSamples,
                                std::uint64_t seed = std::mt19937_64::default_seed);

  /**
   * @brief Constructor of OperationQuadratureMCAdvanced, specifying dimensions
   * and the number of samples to use.
   *
   * @param dimensions dimensionality of this problem
   * @param numberOfSamples Number of Monte Carlo samples
   * @param seed Custom seed (defaults to default seed of mt19937_64)
   */
  OperationQuadratureMCAdvanced(size_t dimensions, size_t numberOfSamples,
                                std::uint64_t seed = std::mt19937_64::default_seed);

  /**
   * Descructor
   */
  virtual ~OperationQuadratureMCAdvanced();

  /**
   * @brief Quadrature using advanced MC in @f$\Omega=[0,1]^d@f$.
   *
   * @param alpha Coefficient vector for current grid
   */
  virtual double doQuadrature(sgpp::base::DataVector& alpha);

  /**
   * @brief Quadrature of an arbitrary function using
   * advanced MC in @f$\Omega=[0,1]^d@f$.
   *
   * @param func The function to integrate
   * @param clientdata Optional data to pass to FUNC
   */
  double doQuadratureFunc(FUNC func, void* clientdata);

  /**
   * @brief Quadrature of the @f$L^2@f$-norm of the error,
   * @f$ ||f(x)-u(x)||_{L^2} @f$, between a given function and the
   * current sparse grid function using
   * advanced MC in @f$\Omega=[0,1]^d@f$.
   *
   * @param func The function @f$f(x)@f$
   * @param clientdata Optional data to pass to FUNC
   * @param alpha Coefficient vector for current grid
   */
  double doQuadratureL2Error(FUNC func, void* clientdata, sgpp::base::DataVector& alpha);

  /**
   * @brief Initialize SampleGenerator for NaiveMC
   */
  void useNaiveMonteCarlo();

  /**
   * @brief Initialize SampleGenerator for StratifiedMC
   *
   * @param n Array of dimension proerties
   */
  void useStratifiedMonteCarlo(std::vector<size_t>& n);

  /**
   * @brief Initialize SampleGenerator for LatinHypercubeMC
   */
  void useLatinHypercubeMonteCarlo();

  /**
   * @brief Initialize SampleGenerator for HaltonSequenceGenerator
   */
  void useQuasiMonteCarloWithHaltonSequences();

  /**
   * @brief Initialize SampleGenerator for SobolSequenceGenerator
   */
  void useQuasiMonteCarloWithSobolSequences();

  /**
   * @brief Initialize SampleGenerator for ScrambledSobolSequenceGenerator
   */
  void useQuasiMonteCarloWithScrambledSobolSequences();

  /**
   * @brief Method returns the total number of samples which can be generated
   * according to the sample generator settings (dimensions and subdivision into strata)
   *
   * @return size_t dimension of samples
   */
  size_t getDimensions();

 protected:
  // Pointer to the grid object
  sgpp::base::Grid* grid;
  // Number of MC samples
  size_t numberOfSamples;
  // number of dimensions (same as in Grid, if given)
  size_t dimensions;

  // seed for the sample generator
  std::uint64_t seed;

  // SampleGenerator Instance
  sgpp::quadrature::SampleGenerator* myGenerator;
};

}  // namespace quadrature
}  // namespace sgpp

#endif /* OPERATIONQUADRATUREMCADVANCED_HPP */

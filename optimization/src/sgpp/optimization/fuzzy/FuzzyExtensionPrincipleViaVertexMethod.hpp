// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLEVIAVERTEXMETHOD_HPP
#define SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLEVIAVERTEXMETHOD_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyExtensionPrincipleViaOptimization.hpp>

#include <vector>

namespace sgpp {
namespace optimization {

/**
 * Zadeh's fuzzy extension principle by the vertex method, where the optimization
 * problems are solved by simply taking the best corners of the confidence intervals.
 */
class FuzzyExtensionPrincipleViaVertexMethod : public FuzzyExtensionPrincipleViaOptimization {
 public:
  /**
   * Constructor.
   *
   * @param f                       function through which to propagate the uncertainties
   * @param numberOfAlphaSegments   number of \f$\alpha\f$ segments
   */
  explicit FuzzyExtensionPrincipleViaVertexMethod(
      const ScalarFunction& f,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS);

  /**
   * Copy constructor.
   *
   * @param other   other fuzzy extension principle
   */
  FuzzyExtensionPrincipleViaVertexMethod(const FuzzyExtensionPrincipleViaVertexMethod& other);

  /**
   * Destructor.
   */
  ~FuzzyExtensionPrincipleViaVertexMethod() override;

  /**
   * @param[out]  clone pointer to cloned object
   */
  void clone(std::unique_ptr<FuzzyExtensionPrinciple>& clone) const override;

 protected:
  /// precomputed storage for the powers of two
  std::vector<size_t> powersOfTwo;
  /// temporary vector used in optimizeForSingleAlphaLevel
  base::DataVector xTmp;

  /**
   * Custom preparation method that is called before the parallelized
   * optimizeForSingleAlphaLevel calls.
   */
  void prepareApply() override;

  /**
   * Solve the minimization/maximization problem for a single \f$\alpha\f$ level.
   *
   * @param[in]   j             index of \f$\alpha\f$ level
   * @param[out]  minimumPoint  minimum point
   * @param[out]  minimumValue  minimum function value
   * @param[out]  maximumPoint  maximum point
   * @param[out]  maximumValue  maximum function value
   */
  void optimizeForSingleAlphaLevel(
      size_t j, base::DataVector& minimumPoint, double& minimumValue,
      base::DataVector& maximumPoint, double& maximumValue) override;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLEVIAVERTEXMETHOD_HPP */

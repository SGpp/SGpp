// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLEVIATRANSFORMATION_HPP
#define SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLEVIATRANSFORMATION_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyExtensionPrincipleViaOptimization.hpp>

#include <vector>

namespace sgpp {
namespace optimization {

/**
 * Zadeh's fuzzy extension principle by the transformation method, where the optimization
 * problems are solved by sampling the optimization domains and taking the
 * best points. The transformation method produces inaccurate results if the
 * fuzzy input intervals are not fuzzy numbers, i.e., if the confidence intervals for
 * \f$\alpha = 1\f$ contain more than one point.
 */
class FuzzyExtensionPrincipleViaTransformation : public FuzzyExtensionPrincipleViaOptimization {
 public:
  /**
   * Constructor.
   *
   * @param f                       function through which to propagate the uncertainties
   * @param numberOfAlphaSegments   number of \f$\alpha\f$ segments
   */
  explicit FuzzyExtensionPrincipleViaTransformation(
      const base::ScalarFunction& f,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS);

  /**
   * Copy constructor.
   *
   * @param other   other fuzzy extension principle
   */
  FuzzyExtensionPrincipleViaTransformation(const FuzzyExtensionPrincipleViaTransformation& other);

  /**
   * Destructor.
   */
  ~FuzzyExtensionPrincipleViaTransformation() override;

  /**
   * @param[out]  clone pointer to cloned object
   */
  void clone(std::unique_ptr<FuzzyExtensionPrinciple>& clone) const override;

 protected:
  /// sampling points
  std::vector<std::vector<base::DataVector>> C;
  /// index scaling vector
  std::vector<size_t> gammaSize;
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

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLEVIATRANSFORMATION_HPP */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLEVIAOPTIMIZATION_HPP
#define SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLEVIAOPTIMIZATION_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyExtensionPrinciple.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/optimization/function/scalar/ScaledScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/ScaledScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/ScaledScalarFunctionHessian.hpp>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>

#include <vector>

namespace sgpp {
namespace optimization {

/**
 * Zadeh's fuzzy extension principle by solving optimization problems
 * for each \f$\alpha\f$ level.
 */
class FuzzyExtensionPrincipleViaOptimization : public FuzzyExtensionPrinciple {
 public:
  /**
   * Constructor.
   * By default, MultiStart is used as optimization algorithm.
   *
   * @param f                       function through which to propagate the uncertainties
   * @param numberOfAlphaSegments   number of \f$\alpha\f$ segments
   */
  explicit FuzzyExtensionPrincipleViaOptimization(
      const ScalarFunction& f,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS);

  /**
   * Constructor with custom optimization algorithm.
   *
   * @param optimizer               optimization algorithm and
   *                                function through which to propagate the uncertainties
   * @param numberOfAlphaSegments   number of \f$\alpha\f$ segments
   */
  explicit FuzzyExtensionPrincipleViaOptimization(
      const optimizer::UnconstrainedOptimizer& optimizer,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS);

  /**
   * Copy constructor.
   *
   * @param other   other fuzzy extension principle
   */
  FuzzyExtensionPrincipleViaOptimization(const FuzzyExtensionPrincipleViaOptimization& other);

  /**
   * Destructor.
   */
  ~FuzzyExtensionPrincipleViaOptimization() override;

  /**
   * @param[out]  clone pointer to cloned object
   */
  void clone(std::unique_ptr<FuzzyExtensionPrinciple>& clone) const override;

 protected:
  /// default optimization algorithm
  optimizer::MultiStart defaultOptimizer;
  /// optimization algorithm
  std::unique_ptr<optimizer::UnconstrainedOptimizer> optimizer;
  /// objective function gradient
  std::unique_ptr<ScalarFunctionGradient> fGradient;
  /// objective function Hessian
  std::unique_ptr<ScalarFunctionHessian> fHessian;
  /// scaled objective function (confidence interval to unit hyper-cube)
  std::unique_ptr<ScalarFunction> fScaled;
  /// scaled objective gradient (confidence interval to unit hyper-cube)
  std::unique_ptr<ScalarFunctionGradient> fGradientScaled;
  /// scaled objective Hessian (confidence interval to unit hyper-cube)
  std::unique_ptr<ScalarFunctionHessian> fHessianScaled;

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

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLEVIAOPTIMIZATION_HPP */

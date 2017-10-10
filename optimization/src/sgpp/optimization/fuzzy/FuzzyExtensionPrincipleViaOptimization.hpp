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
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>

#include <vector>

namespace sgpp {
namespace optimization {

class FuzzyExtensionPrincipleViaOptimization : public FuzzyExtensionPrinciple {
 public:
  static const size_t DEFAULT_NUMBER_OF_ALPHA_SEGMENTS = 10;

  explicit FuzzyExtensionPrincipleViaOptimization(
      const ScalarFunction& f,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS);

  explicit FuzzyExtensionPrincipleViaOptimization(
      const optimizer::UnconstrainedOptimizer& optimizer,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS);

  FuzzyExtensionPrincipleViaOptimization(const FuzzyExtensionPrincipleViaOptimization& other);

  ~FuzzyExtensionPrincipleViaOptimization() override;

  FuzzyInterval* apply(const std::vector<const FuzzyInterval*>& xFuzzy) const override;

  size_t getNumberOfAlphaSegments() const;
  void setNumberOfAlphaSegments(size_t numberOfAlphaSegments);

 protected:
  optimizer::MultiStart defaultOptimizer;
  std::unique_ptr<optimizer::UnconstrainedOptimizer> optimizer;
  std::unique_ptr<ScalarFunctionGradient> fGradient;
  std::unique_ptr<ScalarFunctionHessian> fHessian;
  size_t m;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLEVIAOPTIMIZATION_HPP */

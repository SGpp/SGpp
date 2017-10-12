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

  FuzzyInterval* apply(const std::vector<const FuzzyInterval*>& xFuzzy) override;

  size_t getNumberOfAlphaSegments() const;
  void setNumberOfAlphaSegments(size_t numberOfAlphaSegments);

  const base::DataVector& getAlphaLevels() const;
  const std::vector<base::DataVector>& getOptimizationDomainsLowerBounds() const;
  const std::vector<base::DataVector>& getOptimizationDomainsUpperBounds() const;
  const std::vector<base::DataVector>& getMinimumPoints() const;
  const std::vector<base::DataVector>& getMaximumPoints() const;

 protected:
  optimizer::MultiStart defaultOptimizer;
  std::unique_ptr<optimizer::UnconstrainedOptimizer> optimizer;
  std::unique_ptr<ScalarFunctionGradient> fGradient;
  std::unique_ptr<ScalarFunctionHessian> fHessian;
  size_t m;
  base::DataVector alphaLevels;
  std::vector<base::DataVector> optimizationDomainsUpperBounds;
  std::vector<base::DataVector> optimizationDomainsLowerBounds;
  std::vector<base::DataVector> minimumPoints;
  std::vector<base::DataVector> maximumPoints;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLEVIAOPTIMIZATION_HPP */

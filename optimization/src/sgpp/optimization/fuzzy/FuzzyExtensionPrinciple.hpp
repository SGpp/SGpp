// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLE_HPP
#define SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>
#include <sgpp/optimization/fuzzy/FuzzyInterval.hpp>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>

#include <vector>

namespace sgpp {
namespace optimization {

class FuzzyExtensionPrinciple {
 public:
  static const size_t DEFAULT_NUMBER_OF_ALPHA_SEGMENTS = 10;

  FuzzyExtensionPrinciple(
      const ScalarFunction& f,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS);

  FuzzyExtensionPrinciple(
      const optimizer::UnconstrainedOptimizer& optimizer,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS);

  FuzzyExtensionPrinciple(const FuzzyExtensionPrinciple& other);

  void apply(const std::vector<const FuzzyInterval*>& x, std::unique_ptr<FuzzyInterval>& y) const;

  size_t getNumberOfAlphaSegments() const;
  void setNumberOfAlphaSegments(size_t numberOfAlphaSegments);

 protected:
  optimizer::MultiStart defaultOptimizer;
  std::unique_ptr<optimizer::UnconstrainedOptimizer> optimizer;
  std::unique_ptr<ScalarFunction> f;
  std::unique_ptr<ScalarFunctionGradient> fGradient;
  std::unique_ptr<ScalarFunctionHessian> fHessian;
  size_t m;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLE_HPP */

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

class FuzzyExtensionPrincipleViaTransformation : public FuzzyExtensionPrincipleViaOptimization {
 public:
  explicit FuzzyExtensionPrincipleViaTransformation(
      const ScalarFunction& f,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS);

  FuzzyExtensionPrincipleViaTransformation(const FuzzyExtensionPrincipleViaTransformation& other);

  ~FuzzyExtensionPrincipleViaTransformation() override;

  void clone(std::unique_ptr<FuzzyExtensionPrinciple>& clone) const override;

 protected:
  std::vector<std::vector<base::DataVector>> C;
  std::vector<size_t> gammaSize;
  base::DataVector xTmp;

  void prepareApply() override;

  void optimizeForSingleAlphaLevel(
      size_t j, base::DataVector& minimumPoint, double& minimumValue,
      base::DataVector& maximumPoint, double& maximumValue) override;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLEVIATRANSFORMATION_HPP */

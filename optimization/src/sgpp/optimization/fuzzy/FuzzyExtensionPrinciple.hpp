// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLE_HPP
#define SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/fuzzy/FuzzyInterval.hpp>

#include <vector>

namespace sgpp {
namespace optimization {

class FuzzyExtensionPrinciple {
 public:
  static const size_t DEFAULT_NUMBER_OF_ALPHA_SEGMENTS = 10;

  explicit FuzzyExtensionPrinciple(
      const ScalarFunction& f,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS);

  FuzzyExtensionPrinciple(const FuzzyExtensionPrinciple& other);

  virtual ~FuzzyExtensionPrinciple();

  FuzzyInterval* apply(const std::vector<const FuzzyInterval*>& xFuzzy);

  size_t getNumberOfAlphaSegments() const;
  void setNumberOfAlphaSegments(size_t numberOfAlphaSegments);

  const base::DataVector& getAlphaLevels() const;
  const std::vector<base::DataVector>& getOptimizationDomainsLowerBounds() const;
  const std::vector<base::DataVector>& getOptimizationDomainsUpperBounds() const;
  const std::vector<base::DataVector>& getMinimumPoints() const;
  const base::DataVector& getMinimumValues() const;
  const std::vector<base::DataVector>& getMaximumPoints() const;
  const base::DataVector& getMaximumValues() const;

  virtual void clone(std::unique_ptr<FuzzyExtensionPrinciple>& clone) const = 0;

 protected:
  std::unique_ptr<ScalarFunction> f;
  size_t m;
  base::DataVector alphaLevels;
  std::vector<base::DataVector> optimizationDomainsLowerBounds;
  std::vector<base::DataVector> optimizationDomainsUpperBounds;
  std::vector<base::DataVector> minimumPoints;
  base::DataVector minimumValues;
  std::vector<base::DataVector> maximumPoints;
  base::DataVector maximumValues;

  virtual void prepareApply();

  virtual void optimizeForSingleAlphaLevel(
      size_t j, base::DataVector& minimumPoint, double& minimumValue,
      base::DataVector& maximumPoint, double& maximumValue) = 0;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLE_HPP */

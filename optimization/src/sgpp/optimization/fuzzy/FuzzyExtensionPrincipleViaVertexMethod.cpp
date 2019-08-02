// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyExtensionPrincipleViaVertexMethod.hpp>
#include <sgpp/optimization/fuzzy/InterpolatedFuzzyInterval.hpp>

#include <algorithm>
#include <limits>
#include <vector>

namespace sgpp {
namespace optimization {

FuzzyExtensionPrincipleViaVertexMethod::FuzzyExtensionPrincipleViaVertexMethod(
    const base::ScalarFunction& f,
    size_t numberOfAlphaSegments) :
      FuzzyExtensionPrincipleViaOptimization(f, numberOfAlphaSegments) {
}

FuzzyExtensionPrincipleViaVertexMethod::FuzzyExtensionPrincipleViaVertexMethod(
    const FuzzyExtensionPrincipleViaVertexMethod& other) :
        FuzzyExtensionPrincipleViaVertexMethod(*other.f, other.m) {
  powersOfTwo = other.powersOfTwo;
  xTmp = other.xTmp;
}

FuzzyExtensionPrincipleViaVertexMethod::~FuzzyExtensionPrincipleViaVertexMethod() {}

void FuzzyExtensionPrincipleViaVertexMethod::prepareApply() {
  const size_t d = f->getNumberOfParameters();
  xTmp.resize(d);
  powersOfTwo.resize(d + 1);
  powersOfTwo[0] = 1;

  for (size_t t = 0; t < d; t++) {
    powersOfTwo[t + 1] = 2 * powersOfTwo[t];
  }
}

void FuzzyExtensionPrincipleViaVertexMethod::optimizeForSingleAlphaLevel(
    size_t j, base::DataVector& minimumPoint, double& minimumValue,
    base::DataVector& maximumPoint, double& maximumValue) {
  const size_t d = f->getNumberOfParameters();
  const base::DataVector& lowerBounds = optimizationDomainsLowerBounds[j];
  const base::DataVector& upperBounds = optimizationDomainsUpperBounds[j];

  // calculate minimum and maximum on all vertices of the interval box
  minimumValue = std::numeric_limits<double>::infinity();
  maximumValue = -std::numeric_limits<double>::infinity();

  for (size_t k = 0; k < powersOfTwo[d]; k++) {
    for (size_t t = 0; t < d; t++) {
      xTmp[t] = ((k & powersOfTwo[t]) ? lowerBounds[t] : upperBounds[t]);
    }

    const double fx = f->eval(xTmp);

    if (fx < minimumValue) {
      minimumValue = fx;
      minimumPoint = xTmp;
    }

    if (fx > maximumValue) {
      maximumValue = fx;
      maximumPoint = xTmp;
    }
  }
}

void FuzzyExtensionPrincipleViaVertexMethod::clone(
    std::unique_ptr<FuzzyExtensionPrinciple>& clone) const {
  clone = std::unique_ptr<FuzzyExtensionPrinciple>(
      new FuzzyExtensionPrincipleViaVertexMethod(*this));
}

}  // namespace optimization
}  // namespace sgpp

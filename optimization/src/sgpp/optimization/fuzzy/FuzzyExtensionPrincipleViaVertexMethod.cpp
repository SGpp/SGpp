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
    const ScalarFunction& f,
    size_t numberOfAlphaSegments) :
      FuzzyExtensionPrinciple(f),
      m(numberOfAlphaSegments) {
}

FuzzyExtensionPrincipleViaVertexMethod::FuzzyExtensionPrincipleViaVertexMethod(
    const FuzzyExtensionPrincipleViaVertexMethod& other) :
        FuzzyExtensionPrincipleViaVertexMethod(*other.f, other.m) {
}

FuzzyExtensionPrincipleViaVertexMethod::~FuzzyExtensionPrincipleViaVertexMethod() {}

FuzzyInterval* FuzzyExtensionPrincipleViaVertexMethod::apply(
    const std::vector<const FuzzyInterval*>& xFuzzy) const {
  const size_t d = f->getNumberOfParameters();

  // lower and upper bounds of the interval box
  base::DataVector lowerBounds(d);
  base::DataVector upperBounds(d);

  // temporary evaluation point
  base::DataVector x(d);

  // result data
  base::DataVector xData(2 * m + 2);
  base::DataVector alphaData(2 * m + 2);

  // save powers of two
  std::vector<size_t> powersOfTwo(d + 1);
  powersOfTwo[0] = 1;

  for (size_t t = 0; t < d; t++) {
    powersOfTwo[t + 1] = 2 * powersOfTwo[t];
  }

  // iterate through alphas
  for (size_t j = 0; j <= m; j++) {
    const double alpha = static_cast<double>(j) / static_cast<double>(m);

    // determine input parameter confidence interval,
    // directly changing the optimization domain in fScaled
    for (size_t t = 0; t < d; t++) {
      lowerBounds[t] = xFuzzy[t]->evaluateConfidenceIntervalLowerBound(alpha);
      upperBounds[t] = xFuzzy[t]->evaluateConfidenceIntervalUpperBound(alpha);
    }

    // calculate minimum and maximum on all vertices of the interval box
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();

    for (size_t k = 0; k < powersOfTwo[d]; k++) {
      for (size_t t = 0; t < d; t++) {
        x[t] = ((k & powersOfTwo[t]) ? lowerBounds[t] : upperBounds[t]);
      }

      const double fx = f->eval(x);
      min = std::min(fx, min);
      max = std::max(fx, max);
    }

    // save minimum and maximum in result
    xData[j] = min;
    alphaData[j] = alpha;

    xData[2*m+1-j] = max;
    alphaData[2*m+1-j] = alpha;
  }

  // interpolate between alpha data points
  return new InterpolatedFuzzyInterval(xData, alphaData);
}

size_t FuzzyExtensionPrincipleViaVertexMethod::getNumberOfAlphaSegments() const {
  return m;
}

void FuzzyExtensionPrincipleViaVertexMethod::setNumberOfAlphaSegments(
    size_t numberOfAlphaSegments) {
  m = numberOfAlphaSegments;
}

}  // namespace optimization
}  // namespace sgpp

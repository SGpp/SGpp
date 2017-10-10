// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyExtensionPrincipleViaTransformation.hpp>
#include <sgpp/optimization/fuzzy/InterpolatedFuzzyInterval.hpp>

#include <algorithm>
#include <limits>
#include <vector>

namespace sgpp {
namespace optimization {

FuzzyExtensionPrincipleViaTransformation::FuzzyExtensionPrincipleViaTransformation(
    const ScalarFunction& f,
    size_t numberOfAlphaSegments) :
      FuzzyExtensionPrinciple(f),
      m(numberOfAlphaSegments) {
}

FuzzyExtensionPrincipleViaTransformation::FuzzyExtensionPrincipleViaTransformation(
    const FuzzyExtensionPrincipleViaTransformation& other) :
    FuzzyExtensionPrincipleViaTransformation(*other.f, other.m) {
}

FuzzyExtensionPrincipleViaTransformation::~FuzzyExtensionPrincipleViaTransformation() {}

FuzzyInterval* FuzzyExtensionPrincipleViaTransformation::apply(
    const std::vector<const FuzzyInterval*>& xFuzzy) const {
  const size_t d = f->getNumberOfParameters();

  // lower and upper bounds of the interval box
  base::DataVector lowerBounds(d);
  base::DataVector upperBounds(d);

  // alphas
  base::DataVector alphas(m + 1);

  // temporary evaluation point
  base::DataVector x(d);

  // result data
  base::DataVector xData(2 * m + 2);
  base::DataVector alphaData(2 * m + 2);

  // transformed points
  std::vector<std::vector<base::DataVector>> C(m + 1);

  for (size_t j = m + 1; j-- > 0;) {
    alphas[j] = static_cast<double>(j) / static_cast<double>(m);
    C[j].resize(m - j + 1);

    for (size_t t = 0; t < d; t++) {
      lowerBounds[t] = xFuzzy[t]->evaluateConfidenceIntervalLowerBound(alphas[j]);
      upperBounds[t] = xFuzzy[t]->evaluateConfidenceIntervalUpperBound(alphas[j]);
    }

    for (size_t l = 0; l < m - j + 1; l++) {
      C[j][l].resize(d);

      for (size_t t = 0; t < d; t++) {
        if (l == 0) {
          C[j][l][t] = lowerBounds[t];
        } else if (l == m - j) {
          C[j][l][t] = upperBounds[t];
        } else {
          C[j][l][t] = (C[j + 1][l - 1][t] + C[j + 1][l][t]) / 2.0;
        }
      }
    }
  }

  // iterate through alphas
  for (size_t j = 0; j <= m; j++) {
    // size of gamma vector per dimension
    // (the (m+1-j)^(d-t-1)-vector containing only C[j][l][t])
    std::vector<size_t> gammaSize(d);
    gammaSize[d - 1] = 1;

    for (size_t t = d - 1; t-- > 0;) {
      gammaSize[t] = gammaSize[t + 1] * (m + 1 - j);
    }

    // number of points to check (= (m+1-j)^d)
    const size_t K = gammaSize[0] * (m + 1 - j);

    // calculate minimum and maximum on all vertices of the interval box
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();

    for (size_t k = 0; k < K; k++) {
      for (size_t t = 0; t < d; t++) {
        const size_t l = (k / gammaSize[t]) % (m + 1 - j);
        x[t] = C[j][l][t];
      }

      const double fx = f->eval(x);
      min = std::min(fx, min);
      max = std::max(fx, max);
    }

    // save minimum and maximum in result
    xData[j] = min;
    alphaData[j] = alphas[j];

    xData[2*m+1-j] = max;
    alphaData[2*m+1-j] = alphas[j];
  }

  // interpolate between alpha data points
  return new InterpolatedFuzzyInterval(xData, alphaData);
}

size_t FuzzyExtensionPrincipleViaTransformation::getNumberOfAlphaSegments() const {
  return m;
}

void FuzzyExtensionPrincipleViaTransformation::setNumberOfAlphaSegments(
    size_t numberOfAlphaSegments) {
  m = numberOfAlphaSegments;
}

}  // namespace optimization
}  // namespace sgpp

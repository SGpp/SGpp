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
    const base::ScalarFunction& f,
    size_t numberOfAlphaSegments) :
      FuzzyExtensionPrincipleViaOptimization(f, numberOfAlphaSegments) {
}

FuzzyExtensionPrincipleViaTransformation::FuzzyExtensionPrincipleViaTransformation(
    const FuzzyExtensionPrincipleViaTransformation& other) :
    FuzzyExtensionPrincipleViaOptimization(other),
    C(other.C), gammaSize(other.gammaSize), xTmp(other.xTmp) {
}

FuzzyExtensionPrincipleViaTransformation::~FuzzyExtensionPrincipleViaTransformation() {}

void FuzzyExtensionPrincipleViaTransformation::prepareApply() {
  const size_t d = f->getNumberOfParameters();
  C.resize(m + 1);
  gammaSize.resize(d);
  xTmp.resize(d);

  for (size_t j = m + 1; j-- > 0;) {
    C[j].resize(m - j + 1);

    const base::DataVector& lowerBounds = optimizationDomainsLowerBounds[j];
    const base::DataVector& upperBounds = optimizationDomainsUpperBounds[j];

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
}

void FuzzyExtensionPrincipleViaTransformation::optimizeForSingleAlphaLevel(
    size_t j, base::DataVector& minimumPoint, double& minimumValue,
    base::DataVector& maximumPoint, double& maximumValue) {
  const size_t d = f->getNumberOfParameters();

  // size of gamma vector per dimension
  // (the (m+1-j)^(d-t-1)-vector containing only C[j][l][t])
  gammaSize[d - 1] = 1;

  for (size_t t = d - 1; t-- > 0;) {
    gammaSize[t] = gammaSize[t + 1] * (m + 1 - j);
  }

  // number of points to check (= (m+1-j)^d)
  const size_t K = gammaSize[0] * (m + 1 - j);

  // calculate minimum and maximum on all vertices of the interval box
  minimumValue = std::numeric_limits<double>::infinity();
  maximumValue = -std::numeric_limits<double>::infinity();

  for (size_t k = 0; k < K; k++) {
    for (size_t t = 0; t < d; t++) {
      const size_t l = (k / gammaSize[t]) % (m + 1 - j);
      xTmp[t] = C[j][l][t];
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

void FuzzyExtensionPrincipleViaTransformation::clone(
    std::unique_ptr<FuzzyExtensionPrinciple>& clone) const {
  clone = std::unique_ptr<FuzzyExtensionPrinciple>(
      new FuzzyExtensionPrincipleViaTransformation(*this));
}

}  // namespace optimization
}  // namespace sgpp

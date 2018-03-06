// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/CubicSplineInterpolationEvaluator.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

CubicSplineInterpolationEvaluator::CubicSplineInterpolationEvaluator()
    : evaluationPoint(0.0), basisValues(), xValues(), gridCoefficients() {
  evalConfig.type = CombiEvaluatorTypes::Scalar_CubicSplineInterpolation;
}

CubicSplineInterpolationEvaluator::~CubicSplineInterpolationEvaluator() {}

CubicSplineInterpolationEvaluator::CubicSplineInterpolationEvaluator(
    const CubicSplineInterpolationEvaluator& other)
    : evaluationPoint(other.evaluationPoint),
      basisValues(other.basisValues),
      xValues(other.xValues),
      gridCoefficients(other.gridCoefficients) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_CubicSplineInterpolation;
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>
CubicSplineInterpolationEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>(
      new CubicSplineInterpolationEvaluator(*this));
}

bool CubicSplineInterpolationEvaluator::needsOrderedPoints() { return true; }

bool CubicSplineInterpolationEvaluator::needsParameter() { return true; }

void CubicSplineInterpolationEvaluator::setParameter(const FloatScalarVector& param) {
  evaluationPoint = param.value();
  computeBasisValues();
}

void CubicSplineInterpolationEvaluator::setGridPoints(std::vector<double> const& x) {
  xValues = x;

  if (xValues.size() < 2) {
    if (xValues.size() == 1) {
      basisValues.resize(1);
      basisValues[0] = 1;
    }
    return;
  }

  size_t n = x.size() - 1;

  gridCoefficients = std::vector<std::vector<SplineSet>>(x.size(), std::vector<SplineSet>(n));

  for (size_t j = 0; j < x.size(); ++j) {
    std::vector<double> b(n);
    std::vector<double> d(n);
    std::vector<double> h;

    for (size_t i = 0; i < n; ++i) h.push_back(x[i + 1] - x[i]);

    std::vector<double> alpha;
    for (size_t i = 1; i < n; ++i)
      alpha.push_back(3 * ((i + 1 == j) - (i == j)) / h[i] -
                      3 * ((i == j) - (i - 1 == j)) / h[i - 1]);

    std::vector<double> c(n + 1);
    std::vector<double> l(n + 1);
    std::vector<double> mu(n + 1);
    std::vector<double> z(n + 1);
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for (size_t i = 1; i < n; ++i) {
      l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
      mu[i] = h[i] / l[i];
      z[i] = (alpha[i - 1] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    for (int i = static_cast<int>(n - 1); i >= 0; --i) {
      c[i] = z[i] - mu[i] * c[i + 1];
      b[i] = ((i + 1 == static_cast<int>(j)) - (i == static_cast<int>(j))) / h[i] -
             h[i] * (c[i + 1] + 2 * c[i]) / 3;
      d[i] = (c[i + 1] - c[i]) / 3 / h[i];
    }

    for (size_t i = 0; i < n; ++i) {
      gridCoefficients[j][i].a = (i == j);
      gridCoefficients[j][i].b = b[i];
      gridCoefficients[j][i].c = c[i];
      gridCoefficients[j][i].d = d[i];
    }
  }

  computeBasisValues();
}

void CubicSplineInterpolationEvaluator::computeBasisValues() {
  if (xValues.size() < 2) {
    if (xValues.size() == 1) {
      basisValues.resize(1);
      basisValues[0] = 1;
    }
    return;
  }

  size_t j;
  for (j = 0; j < xValues.size(); j++) {
    if (xValues[j] > evaluationPoint) {
      if (j == 0) j++;
      break;
    }
  }
  j--;

  double dx = evaluationPoint - xValues[j];

  basisValues.resize(gridCoefficients.size());
  for (size_t i = 0; i < gridCoefficients.size(); ++i) {
    basisValues[i] = gridCoefficients[i][j].a + gridCoefficients[i][j].b * dx +
                     gridCoefficients[i][j].c * dx * dx + gridCoefficients[i][j].d * dx * dx * dx;
  }
}

void CubicSplineInterpolationEvaluator::setBasisCoefficientsAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

} /* namespace combigrid */
} /* namespace sgpp */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/L2LejaPointDistribution.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <algorithm>
#include <limits>

namespace sgpp {
namespace combigrid {

void sgpp::combigrid::L2LejaPointDistribution::addPoint(double point) {
  points.push_back(point);
  for (size_t i = 1; i < sortedPoints.size(); ++i) {
    if (sortedPoints[i] > point) {
      sortedPoints.insert(sortedPoints.begin() + i, point);
      return;
    }
  }
}

double L2LejaPointDistribution::evaluate_denominator(size_t& argmaxIndex,
                                                     GaussLegendreQuadrature& quadRule,
                                                     double tol) {
  // ----------------------------------------------------------------
  // evaluate denominator
  double maxLogIntegral = -std::numeric_limits<double>::max();
  argmaxIndex = 0;

  size_t numGaussPoints = points.size() + 1 + numAdditionalPoints;
  base::DataVector roots(numGaussPoints);
  base::DataVector weights(numGaussPoints);

  for (size_t i = 0; i < sortedPoints.size() - 1; ++i) {
    // --------------------------------------------------------------
    // do quadrature
    // performing Gauss-Legendre integration
    numGaussPoints = points.size() + 1 + numAdditionalPoints;

    auto logfunc = [this](double x_unit, double x_prob) {
      double ans = 0.0;
      for (size_t k = 0; k < points.size(); ++k) {
        double diff = x_prob - points[k];
        ans += std::log(diff * diff);
      }
      return ans;
    };

    double a = sortedPoints[i], b = sortedPoints[i + 1];
    double logIntegral = quadRule.evaluate_sumexp_iteratively(
        logfunc, weightFunction, a, b, numGaussPoints, numAdditionalPoints, tol);
    // --------------------------------------------------------------
    if (maxLogIntegral < logIntegral) {
      maxLogIntegral = logIntegral;
      argmaxIndex = i;
    }
  }

  return maxLogIntegral;
}

double L2LejaPointDistribution::evaluate_numerator(size_t argmaxIndex,
                                                   GaussLegendreQuadrature& quadRule, double tol) {
  // performing Gauss-Legendre integration
  size_t numGaussPoints = points.size() + 1 + numAdditionalPoints;

  auto logfunc = [this](double x_unit, double x_prob) {
    double ans = std::log(x_prob);
    for (size_t k = 0; k < points.size(); ++k) {
      double diff = x_prob - points[k];
      ans += std::log(diff * diff);
    }
    return ans;
  };

  double a = sortedPoints[argmaxIndex], b = sortedPoints[argmaxIndex + 1];
  return quadRule.evaluate_sumexp_iteratively(logfunc, weightFunction, a, b, numGaussPoints,
                                              numAdditionalPoints, tol);
}

void L2LejaPointDistribution::computeNextPoint() {
  // initialize Gauss quadrature
  GaussLegendreQuadrature quadRule(points.size() + 1 + numAdditionalPoints);

  // ----------------------------------------------------------------
  // evaluate denominator
  size_t argmaxIndex = 0;
  double logDenominator = evaluate_denominator(argmaxIndex, quadRule);

  // evaluate numerator
  double logNumerator = evaluate_numerator(argmaxIndex, quadRule);
  // --------------------------------------------------------------
  // compute next L2 Leja point
  double x = std::exp(logNumerator - logDenominator);

  if (x >= 1.0 || x < 1e-15) {
    throw base::algorithm_exception(
        "L2LejaPointDistribution::computeNextPoint - no more L2-Leja points available due to "
        "the limited order of numerical quadrature.");
  }

  addPoint(x);
}

L2LejaPointDistribution::L2LejaPointDistribution()
    : points{0.5},
      sortedPoints{0.0, 0.5, 1.0},
      weightFunction(constantFunction<double>(static_cast<double>(1.0))),
      numAdditionalPoints(0) {}

L2LejaPointDistribution::L2LejaPointDistribution(SingleFunction weightFunction,
                                                 size_t numAdditionalPoints)
    : points(),
      sortedPoints{0.0, 1.0},
      weightFunction(weightFunction),
      numAdditionalPoints(numAdditionalPoints) {}

L2LejaPointDistribution::~L2LejaPointDistribution() {}

double L2LejaPointDistribution::compute(size_t numPoints, size_t j) {
  while (points.size() <= j) {
    computeNextPoint();
  }
  return points[j];
}

} /* namespace combigrid */
} /* namespace sgpp */

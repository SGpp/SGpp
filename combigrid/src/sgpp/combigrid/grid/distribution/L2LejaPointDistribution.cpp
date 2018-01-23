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

void L2LejaPointDistribution::computeNextPoint() {
  auto evaluationFunc = [this](double x_unit, double x_trans) {
    long double prod = weightFunction(x_trans);

    for (size_t i = 0; i < points.size(); ++i) {
      prod *= (x_trans - points[i]);
    }

    return prod * prod;
  };

  size_t numQuadPoints = points.size() + 1 + numAdditionalPoints;
  GaussLegendreQuadrature quad(numQuadPoints);

  long double maxIntegral = -1.0;
  size_t argmaxIndex = 0;

  for (size_t i = 0; i < sortedPoints.size() - 1; ++i) {
    long double integral = quad.evaluate_long(evaluationFunc, sortedPoints[i], sortedPoints[i + 1]);

    if (integral > maxIntegral) {
      maxIntegral = integral;
      argmaxIndex = i;
    }
  }

  if (maxIntegral < std::numeric_limits<long double>::min()) {
    throw sgpp::base::algorithm_exception(
        "L2LejaPointDistribution::compute - maximum number of L2-Leja points reached");
  }

  auto funcTimesX = [evaluationFunc](double x_unit, double x_trans) {
    return x_trans * evaluationFunc(x_unit, x_trans);
  };

  long double secondIntegral =
      quad.evaluate_long(funcTimesX, sortedPoints[argmaxIndex], sortedPoints[argmaxIndex + 1]);

  double x = static_cast<double>(secondIntegral / maxIntegral);

  addPoint(x);
}

void L2LejaPointDistribution::computeNextPointLogSumExp() {
  size_t numQuadPoints = points.size() + 1 + numAdditionalPoints;
  auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();
  sgpp::base::DataVector roots, weights;
  quadRule.getLevelPointsAndWeightsNormalized(numQuadPoints, roots, weights);
  GaussLegendreQuadrature quad(numQuadPoints);

  // ----------------------------------------------------------------
  // evaluate denominator
  double maxLogIntegral = -std::numeric_limits<double>::max();
  size_t argmaxIndex = 0;

  sgpp::base::DataVector psi(numQuadPoints);
  for (size_t i = 0; i < sortedPoints.size() - 1; ++i) {
    // --------------------------------------------------------------
    // do quadrature
    double a = sortedPoints[i], b = sortedPoints[i + 1];
    double width = b - a;
    double logIntegral = 0.0;
    double max_psi = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < roots.getSize(); ++i) {
      double x_trans = a + width * roots[i];

      psi[i] = std::log(weights[i]) + std::log(weightFunction(x_trans) * width);
      for (size_t j = 0; j < points.size(); ++j) {
        double diff = x_trans - points[j];
        psi[i] += std::log(diff * diff);
      }

      if (max_psi < psi[i]) {
        max_psi = psi[i];
      }
    }

    double inner_sum = 0.0;
    for (size_t j = 0; j < psi.size(); ++j) {
      inner_sum += std::exp(psi[j] - max_psi);
    }

    logIntegral = max_psi + std::log(inner_sum);
    // --------------------------------------------------------------
    if (maxLogIntegral < logIntegral) {
      maxLogIntegral = logIntegral;
      argmaxIndex = i;
    }
  }

  // ----------------------------------------------------------------
  // compute nominator
  double a = sortedPoints[argmaxIndex], b = sortedPoints[argmaxIndex + 1];
  double width = b - a;
  double max_psi = -std::numeric_limits<double>::max();
  for (size_t i = 0; i < roots.getSize(); ++i) {
    double x_trans = a + width * roots[i];
    psi[i] = std::log(weights[i]) + std::log(weightFunction(x_trans) * width) + std::log(x_trans);

    for (size_t j = 0; j < points.size(); ++j) {
      double diff = x_trans - points[j];
      psi[i] += std::log(diff * diff);
    }

    if (max_psi < psi[i]) {
      max_psi = psi[i];
    }
  }

  double inner_sum = 0.0;
  for (size_t j = 0; j < psi.size(); ++j) {
    inner_sum += std::exp(psi[j] - max_psi);
  }
  double secondLogIntegral = max_psi + std::log(inner_sum);

  // --------------------------------------------------------------
  // compute next L2 Leja point
  double x = std::exp(secondLogIntegral - maxLogIntegral);

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
    computeNextPointLogSumExp();
    //    computeNextPoint();
  }
  return points[j];
}

} /* namespace combigrid */
} /* namespace sgpp */

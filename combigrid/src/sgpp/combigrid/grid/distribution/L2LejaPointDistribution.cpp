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

void L2LejaPointDistribution::evaluate_denominator(size_t& argmaxIndex, double& logIntegral,
                                                   double tol) {
  auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();

  // ----------------------------------------------------------------
  // evaluate denominator
  double maxLogIntegral = -std::numeric_limits<double>::max();
  argmaxIndex = 0;

  for (size_t i = 0; i < sortedPoints.size() - 1; ++i) {
    // --------------------------------------------------------------
    // do quadrature
    double logIntegral_old = 0.0;
    double a = sortedPoints[i], b = sortedPoints[i + 1];
    double width = b - a;

    // performing Gauss-Legendre integration
    size_t numGaussPoints = points.size() + 1 + numAdditionalPoints;
    double err = 1e14;
    size_t iteration = 0;
    base::DataVector roots(numGaussPoints);
    base::DataVector weights(numGaussPoints);

    sgpp::base::DataVector psi(numGaussPoints);

    while (err > tol && numGaussPoints < quadRule.getMaxSupportedLevel()) {
      quadRule.getLevelPointsAndWeightsNormalized(numGaussPoints, roots, weights);
      psi.resize(numGaussPoints);

      double max_psi = -std::numeric_limits<double>::max();
      for (size_t j = 0; j < roots.getSize(); ++j) {
        double x_trans = a + width * roots[j];

        psi[j] = std::log(weights[j]) + std::log(weightFunction(x_trans) * width);
        for (size_t k = 0; k < points.size(); ++k) {
          double diff = x_trans - points[k];
          psi[j] += std::log(diff * diff);
        }

        if (max_psi < psi[j]) {
          max_psi = psi[j];
        }
      }

      double inner_sum = 0.0;
      for (size_t j = 0; j < psi.size(); ++j) {
        inner_sum += std::exp(psi[j] - max_psi);
      }

      logIntegral = max_psi + std::log(inner_sum);

      if (iteration > 0) {
        err = std::fabs(std::exp(logIntegral_old) - std::exp(logIntegral));
      }
      logIntegral_old = logIntegral;
      numGaussPoints += numAdditionalPoints;
      iteration += 1;

      if (numAdditionalPoints == 0) {
        break;
      }
    }

    // --------------------------------------------------------------
    if (maxLogIntegral < logIntegral) {
      maxLogIntegral = logIntegral;
      argmaxIndex = i;
    }
  }

  logIntegral = maxLogIntegral;
}

void L2LejaPointDistribution::evaluate_numerator(size_t argmaxIndex, double& logIntegral,
                                                 double tol) {
  auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();

  // performing Gauss-Legendre integration
  size_t numGaussPoints = points.size() + 1 + numAdditionalPoints;
  double err = 1e14;
  size_t iteration = 0;
  base::DataVector roots(numGaussPoints);
  base::DataVector weights(numGaussPoints);

  sgpp::base::DataVector psi(numGaussPoints);
  double logIntegral_old = 0.0;
  while (err > tol && numGaussPoints < quadRule.getMaxSupportedLevel()) {
    // --------------------------------------------------------------
    // do quadrature
    quadRule.getLevelPointsAndWeightsNormalized(numGaussPoints, roots, weights);
    psi.resize(numGaussPoints);

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
    logIntegral = max_psi + std::log(inner_sum);
    // --------------------------------------------------------------
    if (iteration > 0) {
      err = std::fabs(std::exp(logIntegral_old) - std::exp(logIntegral));
    }
    logIntegral_old = logIntegral;
    numGaussPoints += numAdditionalPoints;
    iteration += 1;

    if (numAdditionalPoints == 0) {
      break;
    }
  }
}

void L2LejaPointDistribution::computeNextPoint() {
  // ----------------------------------------------------------------
  // evaluate denominator
  double firstLogIntegral = 0.0;
  size_t argmaxIndex = 0;
  evaluate_denominator(argmaxIndex, firstLogIntegral);

  // evaluate numerator
  double secondLogIntegral = 0.0;
  evaluate_numerator(argmaxIndex, secondLogIntegral);
  // --------------------------------------------------------------
  // compute next L2 Leja point
  double x = std::exp(secondLogIntegral - firstLogIntegral);

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

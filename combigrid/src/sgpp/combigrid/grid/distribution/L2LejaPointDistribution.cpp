// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/L2LejaPointDistribution.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>

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
  auto evaluationFunc = [this](double x) {
    double prod = weightFunction(x);

    for (size_t i = 0; i < points.size(); ++i) {
      prod *= (x - points[i]);
    }

    return prod * prod;
  };

  size_t numQuadPoints = points.size() + 1 + numAdditionalPoints;
  GaussLegendreQuadrature quad(numQuadPoints);

  double maxIntegral = -1.0;
  size_t argmaxIndex = 0;

  for (size_t i = 0; i < sortedPoints.size() - 1; ++i) {
    double integral = quad.evaluate(evaluationFunc, sortedPoints[i], sortedPoints[i + 1]);

    if (integral > maxIntegral) {
      maxIntegral = integral;
      argmaxIndex = i;
    }
  }

  auto funcTimesX = [evaluationFunc](double x) { return x * evaluationFunc(x); };

  double secondIntegral =
      quad.evaluate(funcTimesX, sortedPoints[argmaxIndex], sortedPoints[argmaxIndex + 1]);

  addPoint(secondIntegral / maxIntegral);
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

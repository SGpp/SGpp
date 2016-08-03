// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <vector>

namespace sgpp {
namespace combigrid {

LinearInterpolationEvaluator::LinearInterpolationEvaluator()
    : evaluationPoint(0.0), basisCoefficients(), xValues() {}

LinearInterpolationEvaluator::~LinearInterpolationEvaluator() {}

LinearInterpolationEvaluator::LinearInterpolationEvaluator(
    const LinearInterpolationEvaluator& other)
    : evaluationPoint(other.evaluationPoint),
      basisCoefficients(other.basisCoefficients),
      xValues(other.xValues) {}

void LinearInterpolationEvaluator::computeBasisCoefficients() {
  size_t numPoints = xValues.size();

  basisCoefficients = std::vector<FloatScalarVector>(
      numPoints, FloatScalarVector(0.0));  // TODO(holzmudd): could be optimized

  if (numPoints == 0) {
    return;
  }

  if (evaluationPoint <= xValues[0]) {
    basisCoefficients[0] = FloatScalarVector::one();
    return;
  }

  // TODO(holzmudd): could be optimized by using binary search
  for (size_t i = 1; i < numPoints; ++i) {
    if (evaluationPoint <= xValues[i]) {
      double x0 = xValues[i - 1];
      double x1 = xValues[i];

      double weight = (x1 - evaluationPoint) / (x1 - x0);
      basisCoefficients[i - 1] = FloatScalarVector(weight);
      basisCoefficients[i] = FloatScalarVector(1 - weight);

      return;
    }
  }

  // if we did not return in the loop, then evaluationPoint > all xValues...
  basisCoefficients[numPoints - 1] = FloatScalarVector::one();
}

void LinearInterpolationEvaluator::setGridPoints(const std::vector<double>& newXValues) {
  xValues = newXValues;
  computeBasisCoefficients();
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >
LinearInterpolationEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >(
      new LinearInterpolationEvaluator(*this));
}

bool LinearInterpolationEvaluator::needsOrderedPoints() { return true; }

bool LinearInterpolationEvaluator::needsParameter() { return true; }

void LinearInterpolationEvaluator::setParameter(const FloatScalarVector& param) {
  evaluationPoint = param.value();
  computeBasisCoefficients();
}

} /* namespace combigrid */
} /* namespace sgpp*/

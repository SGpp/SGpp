// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <vector>

namespace sgpp {
namespace combigrid {

LinearInterpolationEvaluator::LinearInterpolationEvaluator()
    : evaluationPoint(0.0), basisValues(), xValues() {}

LinearInterpolationEvaluator::~LinearInterpolationEvaluator() {}

LinearInterpolationEvaluator::LinearInterpolationEvaluator(
    const LinearInterpolationEvaluator& other)
    : evaluationPoint(other.evaluationPoint),
      basisValues(other.basisValues),
      xValues(other.xValues) {}

void LinearInterpolationEvaluator::computeBasisValues() {
  size_t numPoints = xValues.size();

  basisValues = std::vector<FloatScalarVector>(numPoints, FloatScalarVector(0.0));

  if (numPoints == 0) {
    return;
  }

  if (evaluationPoint <= xValues[0]) {
    if (std::abs(xValues[0]) > 1e-14) {
      basisValues[0] = FloatScalarVector(evaluationPoint / xValues[0]);
    }
    return;
  }

  // TODO(holzmudd): could be optimized by using binary search
  for (size_t i = 1; i < numPoints; ++i) {
    if (evaluationPoint <= xValues[i]) {
      double x0 = xValues[i - 1];
      double x1 = xValues[i];

      double weight = (x1 - evaluationPoint) / (x1 - x0);
      basisValues[i - 1] = FloatScalarVector(weight);
      basisValues[i] = FloatScalarVector(1 - weight);

      return;
    }
  }

  // if we did not return in the loop, then evaluationPoint > all xValues...
  if (std::abs(xValues[numPoints - 1] - 1.0) > 1e-14) {
    basisValues[numPoints - 1] =
        FloatScalarVector((evaluationPoint - 1.0) / (xValues[numPoints - 1] - 1.0));
  }
}

void LinearInterpolationEvaluator::setGridPoints(const std::vector<double>& newXValues) {
  xValues = newXValues;

  // ToDo (rehmemk) Compute B-Spline evaluations at Grid Points. Hand over a flag "SLE" or
  // "FunctionValues" for computeBasisValues()
  computeBasisValues();
}

void LinearInterpolationEvaluator::setFunctionValuesAtGridPoints(
    std::vector<double>& functionValues) {
  // ToDo (rehmemk) Compute coefficients via SLE
  basisCoefficients = functionValues;
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
  computeBasisValues();
}

} /* namespace combigrid */
} /* namespace sgpp*/

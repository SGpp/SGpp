// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/BSplineInterpolationEvaluator.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
#include <sgpp/combigrid/utils/CombigridBSplineBasis.hpp>
#include <vector>
#include <iostream>

namespace sgpp {
namespace combigrid {

BSplineInterpolationEvaluator::BSplineInterpolationEvaluator()
    : evaluationPoint(0.0), basisValues(), xValues(), degree(3) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_BSplineInterpolation;
  evalConfig.degree = 3;
}

BSplineInterpolationEvaluator::BSplineInterpolationEvaluator(size_t degree)
    : evaluationPoint(0.0), basisValues(), xValues(), degree(degree) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_BSplineInterpolation;
  evalConfig.degree = degree;
}

BSplineInterpolationEvaluator::~BSplineInterpolationEvaluator() {}

BSplineInterpolationEvaluator::BSplineInterpolationEvaluator(
    const BSplineInterpolationEvaluator& other)
    : evaluationPoint(other.evaluationPoint),
      basisValues(other.basisValues),
      xValues(other.xValues),
      degree(other.degree) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_BSplineInterpolation;
  evalConfig.degree = other.degree;
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>
BSplineInterpolationEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>(
      new BSplineInterpolationEvaluator(*this));
}

bool BSplineInterpolationEvaluator::needsOrderedPoints() { return true; }

bool BSplineInterpolationEvaluator::needsParameter() { return true; }

void BSplineInterpolationEvaluator::setDegree(size_t const& deg) { degree = deg; }

void BSplineInterpolationEvaluator::setParameter(const FloatScalarVector& param) {
  evaluationPoint = param.value();
  computeBasisValues();
}

void BSplineInterpolationEvaluator::setGridPoints(std::vector<double> const& x) {
  xValues = x;
  computeBasisValues();
}

void BSplineInterpolationEvaluator::computeBasisValues() {
  basisValues.resize(xValues.size(), sgpp::combigrid::FloatScalarVector(0));
  // ToDo (rehmemk) slows down on laptop, test on neon if this is useful
  // #pragma omp parallel for
  for (size_t i = 0; i < xValues.size(); i++) {
    basisValues[i] = expUniformNakBspline(evaluationPoint, degree, i, xValues);
  }
}

void BSplineInterpolationEvaluator::setBasisCoefficientsAtGridPoints(
    std::vector<double>& newBasisCoefficients) {
  basisCoefficients = newBasisCoefficients;
}

} /* namespace combigrid */
} /* namespace sgpp */

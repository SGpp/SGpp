// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <sgpp/combigrid/operation/onedim/ModifiedBSplineInterpolationEvaluator.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
#include <vector>
#include "../../utils/CombigridModifiedNakBSplineBasis.hpp"

namespace sgpp {
namespace combigrid {

ModifiedBSplineInterpolationEvaluator::ModifiedBSplineInterpolationEvaluator()
    : evaluationPoint(0.0), basisValues(), xValues(), degree(3) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_ModifiedBSplineInterpolation;
  evalConfig.degree = 3;
}

ModifiedBSplineInterpolationEvaluator::ModifiedBSplineInterpolationEvaluator(size_t degree)
    : evaluationPoint(0.0), basisValues(), xValues(), degree(degree) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_ModifiedBSplineInterpolation;
  evalConfig.degree = degree;
}

ModifiedBSplineInterpolationEvaluator::~ModifiedBSplineInterpolationEvaluator() {}

ModifiedBSplineInterpolationEvaluator::ModifiedBSplineInterpolationEvaluator(
    const ModifiedBSplineInterpolationEvaluator& other)
    : evaluationPoint(other.evaluationPoint),
      basisValues(other.basisValues),
      xValues(other.xValues),
      degree(other.degree) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_ModifiedBSplineInterpolation;
  evalConfig.degree = other.degree;
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>
ModifiedBSplineInterpolationEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>(
      new ModifiedBSplineInterpolationEvaluator(*this));
}

bool ModifiedBSplineInterpolationEvaluator::needsOrderedPoints() { return true; }

bool ModifiedBSplineInterpolationEvaluator::needsParameter() { return true; }

void ModifiedBSplineInterpolationEvaluator::setDegree(size_t const& deg) { degree = deg; }

void ModifiedBSplineInterpolationEvaluator::setParameter(const FloatScalarVector& param) {
  evaluationPoint = param.value();
  computeBasisValues();
}

void ModifiedBSplineInterpolationEvaluator::setGridPoints(std::vector<double> const& x) {
  xValues = x;
  computeBasisValues();
}

void ModifiedBSplineInterpolationEvaluator::computeBasisValues() {
  basisValues.resize(xValues.size(), sgpp::combigrid::FloatScalarVector(0));
  // ToDo (rehmemk) slows down on laptop, test on neon if this is useful
  // #pragma omp parallel for
  for (size_t i = 0; i < xValues.size(); i++) {
    basisValues[i] = expUniformModifiedNakBspline(evaluationPoint, degree, i, xValues);
  }
}

void ModifiedBSplineInterpolationEvaluator::setBasisCoefficientsAtGridPoints(
    std::vector<double>& newBasisCoefficients) {
  basisCoefficients = newBasisCoefficients;
}

} /* namespace combigrid */
} /* namespace sgpp */

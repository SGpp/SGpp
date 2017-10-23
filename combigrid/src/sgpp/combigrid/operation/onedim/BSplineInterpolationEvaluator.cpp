// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/BSplineInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

BSplineInterpolationEvaluator::BSplineInterpolationEvaluator()
    : evaluationPoint(0.0), basisValues(), xValues(), degree(3) {}

BSplineInterpolationEvaluator::BSplineInterpolationEvaluator(size_t degree)
    : evaluationPoint(0.0), basisValues(), xValues(), degree(degree) {}

BSplineInterpolationEvaluator::~BSplineInterpolationEvaluator() {}

BSplineInterpolationEvaluator::BSplineInterpolationEvaluator(
    const BSplineInterpolationEvaluator& other)
    : evaluationPoint(other.evaluationPoint),
      basisValues(other.basisValues),
      xValues(other.xValues),
      degree(other.degree) {}

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
  /*
   * Am Rand: Spiegle Gitterweiten nach außen
   *
   * nonuniform B-Spline vom Grad n:
   * Knotenfolge xi_k,.,xi_{k+n+1}
   * b^n_{k,xi} = gamma b^{n-1}_{k,xi} + (1-gamma) b^{n-1}_{k+1,xi}
   * gamma(x) = (x-xi_k)/(xi_{k+n}-xi_k)
   *
   *
   */

  basisValues.resize(xValues.size(), sgpp::combigrid::FloatScalarVector(0));

  // Testing:
  //  for (size_t i = 0; i < xValues.size(); i++) {
  //	  //this is f(x) = x
  //    basisValues[i] = evaluationPoint ;
  //  }
  //  return;
  //**************************

  if (xValues.size() == 1) {
    basisValues[0] = 1.0;
    return;
  }
  // Lagrange polynomials for less than 9 points because 9 is the number of gridpoints of a uniform
  // boundary grid of level 3 and this is the first level with enough gridpoints for nak B-Splines
  // Should work for degree 5 as well
  // For degree 7 and higher level 3 with nak is too small to provide enough knots even for one
  // single spline
  else if (xValues.size() < 9) {
    for (size_t i = 0; i < xValues.size(); i++) {
      basisValues[i] = LagrangePolynomial(evaluationPoint, xValues, i);
    }
    return;
  }
  std::vector<double> xi;
  createNakKnots(xValues, degree, xi);
  for (size_t i = 0; i < xValues.size(); i++) {
    basisValues[i] = nonUniformBSpline(evaluationPoint, degree, i, xi);
  }
}

void BSplineInterpolationEvaluator::setFunctionValuesAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

} /* namespace combigrid */
} /* namespace sgpp */

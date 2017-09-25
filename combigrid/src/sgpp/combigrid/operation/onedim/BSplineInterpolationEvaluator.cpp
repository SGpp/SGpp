// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/BSplineInterpolationEvaluator.hpp>

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

BSplineInterpolationEvaluator::BSplineInterpolationEvaluator()
    : evaluationPoint(0.0), basisValues(), xValues(), degree(3) {}

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

  //  if (xValues.size() < 2) {
  //    if (xValues.size() == 1) {
  //      basisValues.resize(1);
  //      basisValues[0] = 1;
  //    }
  //    return;
  //  }

  computeBasisValues();
}

/**
   * @param x     evaluation point
   * @param deg     B-spline degree
   * @param k     index of B-spline in the knot sequence
   * @return      value of non-uniform B-spline
   *              with knots \f$\{\xi_k, ... \xi_{k+p+1}\}\f$
   */
double BSplineInterpolationEvaluator::nonUniformBSpline(double x, size_t deg, size_t k) {
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1
  std::vector<double> xi(degree + 1, 0);
  // this offset is only correct for odd degrees
  size_t offset = (degree + 1) / 2;
  xi.insert(xi.begin() + offset, xValues.begin(), xValues.end());
  for (size_t i = 0; i < offset; i++) {
    xi[offset - i - 1] = -xValues[i + 1] + 2 * xValues[0];
    xi[xValues.size() + offset + i] =
        xValues[xValues.size() - 1] +
        (xValues[xValues.size() - 1] - xValues[xValues.size() - i - 2]);
  }

  //  for (size_t i = 0; i < xi.size(); i++) {
  //    std::cout << xi[i] << " ";
  //  }
  //  std::cout << "\n";

  if (deg == 0) {
    // characteristic function of [xi[k], xi[k+1])
    return (((x >= xi[k]) && (x < xi[k + 1])) ? 1.0 : 0.0);
  } else if ((x < xi[k]) || (x >= xi[k + deg + 1])) {
    // out of support
    return 0.0;
  } else {
    // Cox-de-Boor recursion
    return (x - xi[k]) / (xi[k + deg] - xi[k]) * nonUniformBSpline(x, deg - 1, k) +
           (1.0 - (x - xi[k + 1]) / (xi[k + deg + 1] - xi[k + 1])) *
               nonUniformBSpline(x, deg - 1, k + 1);
  }
}

void BSplineInterpolationEvaluator::computeBasisValues() {
  /*mögliche Belegungen von xValues:
   *
   * 0.5
   * 0 0.5 1
   * 0 0.25 0.5 0.75 1
   * 0 0.125 0.25 0.375 0.5 0.625 0.75 0.875 1
   * ...
   *
   * Im Inneren: nonuniform B-Splines
   * Heuristiken (Andere Ansätze möglich[Mehrfachpunkte, andere Hilfspunkte außen]):
   * Am Rand: Spiegle Gitterweiten nach außen
   * Nur ein Punkt: Setze Giterweite auf 1
   * Ab zwei Punkten -> Randfall
   *
   * nonuniform B-Spline vom Grad n:
   * Knotenfolge xi_k,.,xi_{k+n+1}
   * b^n_{k,xi} = gamma b^{n-1}_{k,xi} + (1-gamma) b^{n-1}_{k+1,xi}
   * gamma(x) = (x-xi_k)/(xi_{k+n}-xi_k)
   *
   *
   */

  //  for (size_t i = 0; i < xValues.size(); i++) {
  //    std::cout << xValues[i] << " ";
  //  }
  //  std::cout << "\n";

  // ToDo(rehmemk) Diese Spezialfälle behandeln
  if (xValues.size() == 1) {
    basisValues.resize(1);
    basisValues[0] = 1;
    return;
  }
  basisValues.resize(xValues.size(), sgpp::combigrid::FloatScalarVector(0));
  //  if (xValues.size() < degree + 1) {
  //    for (size_t i = 0; i < xValues.size(); i++) {
  //      basisValues[i] = 1;
  //    }
  //    return;
  //  }

  // ToDo (rehmemk) Randfälle !
  for (size_t i = 0; i < xValues.size(); i++) {
    //    std::vector<double> xi(degree + 2);
    //    for (size_t k = 0; k < degree + 2; k++) {
    //      xi[k] = xValues[i + k - (degree + 1) / 2];
    //    }
    //    std::cout << i << ": ";
    //    for (size_t l = 0; l < xi.size(); l++) {
    //      std::cout << xi[l] << " ";
    //    }
    //    std::cout << "\n";
    basisValues[i] = nonUniformBSpline(evaluationPoint, degree, i);
  }
}

void BSplineInterpolationEvaluator::setFunctionValuesAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

} /* namespace combigrid */
} /* namespace sgpp */

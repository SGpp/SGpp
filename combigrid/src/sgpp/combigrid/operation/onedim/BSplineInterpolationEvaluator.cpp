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
  // ToDo(rehmemk) this offset is only correct for odd degrees
  size_t offset = (degree + 1) / 2;
  std::vector<double> xi(2 * offset, 0);

  // ToDo(rehmemk) Don't calculate xi every single time!
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

/**
   * @param x     evaluation point
   * @param deg     B-spline degree
   * @param k     index of B-spline in the knot sequence
   * @return      value of non-uniform B-spline
   *              with knots \f$\{\xi_k, ... \xi_{k+p+1}\}\f$
   */
double BSplineInterpolationEvaluator::nonUniformNAKBSpline(double x, size_t deg, size_t k) {
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1
  // ToDo(rehmemk) this offset is only correct for odd degrees
  size_t offset = (degree + 1) / 2;
  std::vector<double> xi(2 * offset + 2, 0);

  // ToDo(rehmemk) Don't calculate xi every single time!
  xi.insert(xi.begin() + offset + 1, xValues.begin(), xValues.end());

  for (size_t i = 0; i < offset + 1; i++) {
    xi[offset - i] = -xValues[i + 1] + 2 * xValues[0];
    xi[xValues.size() + offset + i + 1] =
        xValues[xValues.size() - 1] +
        (xValues[xValues.size() - 1] - xValues[xValues.size() - i - 2]);
  }

  xi.erase(xi.begin() + offset + 2);
  xi.erase(xi.end() - offset - 3);

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

/**
   * @param x     evaluation point
   * @param k     index in the knot sequence
   * @return      value of Lagrange Polynomial
   */
double BSplineInterpolationEvaluator::LagrangePolynomial(double x, size_t k) {
  double res = 1.0;
  for (size_t m = 0; m < xValues.size(); m++) {
    if (k != m) {
      res *= (x - xValues[m]) / (xValues[k] - xValues[m]);
    }
  }
  return res;
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
  if (xValues.size() == 1) {
    basisValues[0] = 1.0;
    return;
  }
  // 9 weil das der Menge an Gitterpunkten eines UniformBoundaryGrid Level 3 entspricht.
  // Todo (rehmemk) Überlegen, wie sich das mit den anderen Gittertypen verhält und nak für kleine
  // Level entsprechend anpassen. Lagrange Polynome!
  else if (xValues.size() < 9) {
    for (size_t i = 0; i < xValues.size(); i++) {
      basisValues[i] = LagrangePolynomial(evaluationPoint, i);
    }
    return;
  }

  for (size_t i = 0; i < xValues.size(); i++) {
    basisValues[i] = nonUniformNAKBSpline(evaluationPoint, degree, i);
  }
}

void BSplineInterpolationEvaluator::setFunctionValuesAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

} /* namespace combigrid */
} /* namespace sgpp */

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
  computeBasisValues();
}

/**
   * @param x     evaluation point
   * @param deg     B-spline degree
   * @param k     index of B-spline in the knot sequence
   * @param xi    vector containing the B-Splines knots
   * @return      value of non-uniform B-spline
   *              with knots \f$\{\xi_k, ... \xi_{k+p+1}\}\f$
   */
double BSplineInterpolationEvaluator::nonUniformBSpline(double x, size_t deg, size_t k,
                                                        std::vector<double> const& xi) {
  if (deg == 0) {
    // characteristic function of [xi[k], xi[k+1])
    return (((x >= xi[k]) && (x < xi[k + 1])) ? 1.0 : 0.0);
  } else if ((x < xi[k]) || (x >= xi[k + deg + 1])) {
    // out of support
    return 0.0;
  } else {
    // Cox-de-Boor recursion
    return (x - xi[k]) / (xi[k + deg] - xi[k]) * nonUniformBSpline(x, deg - 1, k, xi) +
           (1.0 - (x - xi[k + 1]) / (xi[k + deg + 1] - xi[k + 1])) *
               nonUniformBSpline(x, deg - 1, k + 1, xi);
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

/**
   * @param xi vector containing the knots with which the Bsplines are created. This is the most
 * simple case. xi = x inside [0,1] and at the left and right end the necessary amount of inner
 * points are mirrored to the outside
   */
void BSplineInterpolationEvaluator::createKnots(std::vector<double>& xi) {
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1
  // ToDo(rehmemk) this offset is only correct for odd degrees
  size_t offset = (degree + 1) / 2;
  xi.resize(2 * offset, 0);
<<<<<<< HEAD
  //  for (size_t i = 0; i < xi.size(); i++) {
  //    std::cout << xi[i] << " ";
  //  }
  //  std::cout << "\n";

  xi.insert(xi.begin() + offset, xValues.begin(), xValues.end());
  //  for (size_t i = 0; i < xi.size(); i++) {
  //    std::cout << xi[i] << " ";
  //  }
  //  std::cout << "\n";
=======
  xi.insert(xi.begin() + offset, xValues.begin(), xValues.end());

>>>>>>> origin/newCombigridModule
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
}

/**
   * @param xi vector containing the knots with which the Bsplines are created. For dealing with
   * the boundaries at 0 and 1 not a knot knots are used. In the case of degree 3 this means that
   * the knot directly to the right/left of 0/1 are removed.
   */
void BSplineInterpolationEvaluator::createdeg3NakKnots(std::vector<double>& xi) {
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1
<<<<<<< HEAD
  if (degree == 1) {
    createKnots(xi);
    //    std::cout << "The not a knot concept is not applicable for B-Splines of degree 1! "
    //                 "Regular knots will be used for this calculation."
    //              << std::endl;
    return;
  }
=======
>>>>>>> origin/newCombigridModule

  // ToDo(rehmemk) this offset is only correct for odd degrees
  size_t offset = (degree + 1) / 2;
  xi.resize(2 * offset + 2, 0);

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
}

/**
   * @param xi vector containing the knots with which the Bsplines are created. For dealing with
   * the boundaries at 0 and 1 not a knot knots are used. In the case of degree 5 this means that
   * the two knots directly to the right/left of 0/1 are removed.
   */
void BSplineInterpolationEvaluator::createdeg5NakKnots(std::vector<double>& xi) {
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1

  // ToDo(rehmemk) this offset is only correct for odd degrees
  size_t offset = (degree + 1) / 2;
  xi.resize(2 * (offset + 2), 0);

  // ToDo(rehmemk) Don't calculate xi every single time!
  xi.insert(xi.begin() + offset + 2, xValues.begin(), xValues.end());

  for (size_t i = 0; i < offset + 2; i++) {
    xi[offset - i + 1] = -xValues[i + 1] + 2 * xValues[0];
    xi[xValues.size() + offset + i + 2] =
        xValues[xValues.size() - 1] +
        (xValues[xValues.size() - 1] - xValues[xValues.size() - i - 2]);
  }

  xi.erase(xi.begin() + offset + 3);
  xi.erase(xi.begin() + offset + 3);
  xi.erase(xi.end() - offset - 4);
  xi.erase(xi.end() - offset - 4);
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
  // Lagrange polynomials for less than 9 points because 9 is the number of gridpoints of a uniform
  // boundary grid of level 3 and this is the first level with enough gridpoints for nak B-Splines
  // (this is a very heuristic motivation but it works so far. Feel free to implement something
<<<<<<< HEAD
  // better)
  // Should work for degree 5 as well. For degree 7 and higher level 3 with nak is too small to
=======
  // better). SHould work for degree 5 as well
  // For degree 7 and higher level 3 with nak is too small to
>>>>>>> origin/newCombigridModule
  // provide enough knots even for one single spline
  else if (xValues.size() < 9) {
    for (size_t i = 0; i < xValues.size(); i++) {
      basisValues[i] = LagrangePolynomial(evaluationPoint, i);
    }
    return;
  }
  std::vector<double> xi;
  // Choose between nak (not a knot) knots or regular knots
<<<<<<< HEAD
  createKnots(xi);
  //  createNakKnots(xi);
=======
  if (degree == 1) {
    createKnots(xi);
  } else if (degree == 3) {
    createdeg3NakKnots(xi);
  } else if (degree == 5) {
    createdeg5NakKnots(xi);
  } else {
    throw std::invalid_argument("only B-Spline degrees 1,3 and 5 supported");
  }
>>>>>>> origin/newCombigridModule
  for (size_t i = 0; i < xValues.size(); i++) {
    basisValues[i] = nonUniformBSpline(evaluationPoint, degree, i, xi);
    //    std::cout << xValues[i] << " ";
  }
  //  std::cout << "\n";
  //  for (size_t i = 0; i < xi.size(); i++) {
  //    std::cout << xi[i] << " ";
  //  }
  //  std::cout << "\n";
}

void BSplineInterpolationEvaluator::setFunctionValuesAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

} /* namespace combigrid */
} /* namespace sgpp */

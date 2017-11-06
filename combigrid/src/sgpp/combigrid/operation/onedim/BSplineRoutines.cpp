// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>
/**
   * @param x     evaluation point
   * @param deg     B-spline degree
   * @param k     index of B-spline in the knot sequence
   * @param xi    vector containing the B-Splines knots
   * @return      value of non-uniform B-spline
   *              with knots \f$\{\xi_k, ... \xi_{k+p+1}\}\f$
   */
double nonUniformBSpline(double const& x, size_t const& deg, size_t const& index,
                         std::vector<double> const& xi) {
  if (deg == 0) {
    // characteristic function of [xi[k], xi[k+1])
    return (((x >= xi[index]) && (x < xi[index + 1])) ? 1.0 : 0.0);
  } else if ((x < xi[index]) || (x >= xi[index + deg + 1])) {
    // out of support
    return 0.0;
  } else {
    // Cox-de-Boor recursion
    return (x - xi[index]) / (xi[index + deg] - xi[index]) *
               nonUniformBSpline(x, deg - 1, index, xi) +
           (1.0 - (x - xi[index + 1]) / (xi[index + deg + 1] - xi[index + 1])) *
               nonUniformBSpline(x, deg - 1, index + 1, xi);
  }
}

/**
   * @param x     evaluation point
   * @param k     index in the knot sequence
   * @return      value of Lagrange Polynomial
   */
double LagrangePolynomial(double const& x, std::vector<double> const& xValues, size_t const& k) {
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
void createKnots(std::vector<double> const& xValues, size_t const& degree,
                 std::vector<double>& xi) {
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1

  // this offset works only for odd degrees. if even degrees shall be supported it must be
  // generalized
  size_t offset = (degree + 1) / 2;
  xi.resize(2 * offset, 0);
  xi.insert(xi.begin() + offset, xValues.begin(), xValues.end());

  for (size_t i = 0; i < offset; i++) {
    xi[offset - i - 1] = -xValues[i + 1] + 2 * xValues[0];
    xi[xValues.size() + offset + i] =
        xValues[xValues.size() - 1] +
        (xValues[xValues.size() - 1] - xValues[xValues.size() - i - 2]);
  }
}

/**
   * @param xi vector containing the knots with which the Bsplines are created. For dealing with
   * the boundaries at 0 and 1 not a knot knots are used. In the case of degree 3 this means that
   * the knot directly to the right/left of 0/1 are removed.
   */
void createdeg3NakKnots(std::vector<double> const& xValues, size_t const& degree,
                        std::vector<double>& xi) {
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1

  size_t offset = (degree + 1) / 2;
  xi.resize(2 * offset + 2, 0);

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
void createdeg5NakKnots(std::vector<double> const& xValues, size_t const& degree,
                        std::vector<double>& xi) {
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1

  size_t offset = (degree + 1) / 2;
  xi.resize(2 * (offset + 2), 0);

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

void createNakKnots(std::vector<double> const& xValues, size_t const& degree,
                    std::vector<double>& xi) {
  if (degree == 1) {
    createKnots(xValues, degree, xi);
  } else if (degree == 3) {
    createdeg3NakKnots(xValues, degree, xi);
  } else if (degree == 5) {
    createdeg5NakKnots(xValues, degree, xi);
  } else {
    throw std::invalid_argument("BSplineRoutines: only B-Spline degrees 1,3 and 5 supported");
  }
}

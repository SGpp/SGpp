// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>

#include <algorithm>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * log 2 for size_t types
 * @param n size_t number
 * @return log2 of the parameter as size_t
 */
size_t log2(size_t n);

/**
* evaluates a not a knot Bspline on an expUnifromGrid, given by its degree, index and the knot
* sequence it is defined on in x. This routine is much faster than the general nonUniformBSpline.
*@param x       evaluation point
*@param degree     B-spline degree
*@param i       index of B-spline
*@param points    points of the 1D grid
*@return        value of non-uniform B-spline in x
*/
double expUniformNakBspline(double const& x, size_t const& degree, size_t i,
                            std::vector<double> const& points);

/**
 * evaluates a Bspline given by its degree, index and the knot sequence it is defined on in x
   * @param x     evaluation point
   * @param deg     B-spline degree
   * @param index     index of B-spline in the knot sequence
   * @param xi    vector containing the B-Splines knots
   * @return      value of non-uniform B-spline in x
   */
double nonUniformBSpline(double const& x, size_t const& deg, size_t const& index,
                         std::vector<double> const& xi);

/**
 * evaluates the Lagrange polynomial given by its index and the knot sequence it is defined on in x
   * @param x     evaluation point
   * @param xValues
   * @param k     index in the knot sequence
   * @return      value of Lagrange polynomial in x
   */
double LagrangePolynomial(double const& x, std::vector<double> const& xValues, size_t const& k);
}  // namespace combigrid
}  // namespace sgpp

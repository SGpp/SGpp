// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_BSPLINEROUTINES_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_BSPLINEROUTINES_HPP_

#include <sgpp/combigrid/definitions.hpp>

/**
   * @param x     evaluation point
   * @param deg     B-spline degree
   * @param k     index of B-spline in the knot sequence
   * @param xi    vector containing the B-Splines knots
   * @return      value of non-uniform B-spline
   *              with knots \f$\{\xi_k, ... \xi_{k+p+1}\}\f$
   */
double nonUniformBSpline(double const& x, size_t const& deg, size_t const& k,
                         std::vector<double> const& xi);

/**
   * @param x     evaluation point
   * @param k     index in the knot sequence
   * @return      value of Lagrange Polynomial
   */
double LagrangePolynomial(double const& x, std::vector<double> const& xValues, size_t const& k);

/**
   * @param xi vector containing the knots with which the Bsplines are created. This is the most
 * simple case. xi = x inside [0,1] and at the left and right end the necessary amount of inner
 * points are mirrored to the outside
   */
void createKnots(std::vector<double> const& xValues, size_t const& degree, std::vector<double>& xi);

/**
   * @param xi vector containing the knots with which the Bsplines are created. For dealing with
   * the boundaries at 0 and 1 not a knot knots are used. In the case of degree 3 this means that
   * the knot directly to the right/left of 0/1 are removed.
   */
void createdeg3NakKnots(std::vector<double> const& xValues, size_t const& degree,
                        std::vector<double>& xi);

/**
   * @param xi vector containing the knots with which the Bsplines are created. For dealing with
   * the boundaries at 0 and 1 not a knot knots are used. In the case of degree 5 this means that
   * the two knots directly to the right/left of 0/1 are removed.
   */
void createdeg5NakKnots(std::vector<double> const& xValues, size_t const& degree,
                        std::vector<double>& xi);

void createNakKnots(std::vector<double> const& xValues, size_t const& degree,
                    std::vector<double>& xi);
#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_BSPLINEROUTINES_HPP_ */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_BSPLINEROUTINES_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_BSPLINEROUTINES_HPP_

#include <sgpp/combigrid/definitions.hpp>

#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

typedef sgpp::combigrid::GeneralFunction<double, sgpp::base::DataVector const&> MultiFunction;

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

/**
   * @param numDimensions   number of dimensions
   * @param func     		objective function
   * @param degree       	B-spline degree
   * * @return      		Operation for calculating mean(u^2) where u is the B-Spline interpolant
 * of
   * 						func
   * Mean(u^2) = \int u^2(x) f(x) dx = \sum_i \alpha_i \sum_j \alpha_j \int b_i b_j dx
   * Therefore this operation uses the special quadrature for \int b_i b_j and a custom grid
 * function that calculated the interpolation coefficients \alpha and saves the products \alpha_i
 * \alpha_j into a TreeStorage
   */
std::shared_ptr<sgpp::combigrid::CombigridOperation> BsplineMeanSquare(size_t numDimensions,
                                                                       MultiFunction func,
                                                                       size_t degree);
#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_BSPLINEROUTINES_HPP_ */

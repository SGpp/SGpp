// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_BSPLINEROUTINES_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_BSPLINEROUTINES_HPP_

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluationStrategy.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <vector>

/**
* evaluates a not a knot Bspline on an expUnifromGrid, given by its degree, index and the knot
* sequence it is defined on in x. This routine is much faster than the general nonUniformBSpline.
*@param x       evaluation point
*@param degree     B-spline degree
*@param i       index of B-spline
*@param points		points of the 1D grid
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

/**
 * Creates the knot sequence xi needed for the evaluation of B-splines from the evaluation points
 * xValues for B splines of degree 1 by adding the necessary point outisde [0,1] by mirroring at 0
 * and 1.
 * @param xValues grid points inside [0,1]
 * @param xi returns the vector containing the knots with which the Bsplines are created.
 */
std::vector<double> createdeg1Knots(std::vector<double> const& xValues);

/**
 * Creates the knot sequence xi needed for the evaluation of B-splines from the evaluation points
 * xValues for B splines of degree 1 by adding the necessary point outisde [0,1] by mirroring at 0
 * and 1. For dealing with the boundaries at 0 and 1 not a knot knots are used. In the case of
 * degree 3 this means that the knot directly to the right/left of 0/1 are removed.
 * @param xValues grid points inside [0,1]
 * @param xi returns the vector containing the knots with which the Bsplines are created.
   */
std::vector<double> createdeg3NakKnots(std::vector<double> const& xValues);

/**
 * Creates the knot sequence xi needed for the evaluation of B-splines from the evaluation points
 * xValues for B splines of degree 1 by adding the necessary point outisde [0,1] by mirroring at 0
 * and 1. For dealing with the boundaries at 0 and 1 not a knot knots are used. In the case of
 * degree 3 this means that the two knots directly to the right/left of 0/1 are removed.
 * @param xValues grid points inside [0,1]
 * @param xi returns the vector containing the knots with which the Bsplines are created.
   */
std::vector<double> createdeg5NakKnots(std::vector<double> const& xValues);
/**
 * interface for creating the knot sequence xi needed for the evaluation of B-splines from the
 * evaluation points
 * xValues for B splines of degrees i n{1,3,5}
 * @param xValues grid points inside [0,1]
 * @param degree degree of the B spline basis functions
 * @param xi returns the vector containing the knots with which the Bsplines are created.
   */
std::vector<double> createNakKnots(std::vector<double> const& xValues, size_t const& degree);

/**
 * Creates the GridFunction that calculates the coefficients of the B-spline interpolation.
 * The coefficients for each B-Spline are saved in a TreeStorage encoded by a MultiIndex
 * The Grid Functions coefficients are used as well for quadrature.
 *
 * @param func the objective function that shall be interpolated
 * @param grids vector of one dimensional grids
 * @param degree degree of the B spline basis functions
 */
// The coefficients for each B-Spline are saved in a TreeStorage encoded by a MultiIndex
sgpp::combigrid::GridFunction BSplineCoefficientGridFunction(
    sgpp::combigrid::MultiFunction func, sgpp::combigrid::CombiHierarchies::Collection grids,
    size_t degree);

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_BSPLINEROUTINES_HPP_ */

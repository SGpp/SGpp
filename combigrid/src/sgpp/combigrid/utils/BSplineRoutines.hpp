// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/common/GridConversion.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/LevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluationStrategy.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/quadrature/sampling/NaiveSampleGenerator.hpp>

#include <vector>

#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/system/FullSLE.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>

/**
 * Creates the knot sequence xi needed for the evaluation of B-splines from the evaluation points
 * xValues for B splines of degree 1 by adding the necessary point outisde [0,1] by mirroring at 0
 * and 1.
 * @param xValues grid points inside [0,1]
 */
std::vector<double> createdeg1Knots(std::vector<double> const& xValues);

/**
 * Creates the knot sequence xi needed for the evaluation of B-splines from the evaluation points
 * xValues for B splines of degree 1 by adding the necessary point outisde [0,1] by mirroring at 0
 * and 1. For dealing with the boundaries at 0 and 1 not a knot knots are used. In the case of
 * degree 3 this means that the knot directly to the right/left of 0/1 are removed.
 * @param xValues grid points inside [0,1]
   */
std::vector<double> createdeg3NakKnots(std::vector<double> const& xValues);

/**
 * Creates the knot sequence xi needed for the evaluation of B-splines from the evaluation points
 * xValues for B splines of degree 1 by adding the necessary point outisde [0,1] by mirroring at 0
 * and 1. For dealing with the boundaries at 0 and 1 not a knot knots are used. In the case of
 * degree 3 this means that the two knots directly to the right/left of 0/1 are removed.
 * @param xValues grid points inside [0,1]
   */
std::vector<double> createdeg5NakKnots(std::vector<double> const& xValues);
/**
 * interface for creating the knot sequence xi needed for the evaluation of B-splines from the
 * evaluation points
 * xValues for B splines of degrees i n{1,3,5}
 * @param xValues grid points inside [0,1]
 * @param degree degree of the B spline basis functions
   */
std::vector<double> createNakKnots(std::vector<double> const& xValues, size_t const& degree);

/**
 * unique index for level index pair
 *
 * @param level
 * @param index
 * @return
 */
size_t getUniqueIndex(size_t level, size_t index);

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
sgpp::combigrid::GridFunction BSplineTensorCoefficientGridFunction(
    sgpp::combigrid::MultiFunction func, sgpp::combigrid::CombiHierarchies::Collection grids,
    size_t degree);

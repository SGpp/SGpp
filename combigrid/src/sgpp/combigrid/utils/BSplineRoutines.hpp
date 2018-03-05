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
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/quadrature/sampling/NaiveSampleGenerator.hpp>

#include <vector>

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

/**
 * Creates a level structure according to an averaging level manager using variance calculations on
 * each level as norm. This is a very specific case created for the CO2 example. It can (should?) be
 * generalized
 *
 * @param degree        B spline degree
 * @param numDimensions number of dimensions
 * @param func		      the objective function
 * @param levelManager  level manager
 * @return a combigrid operation calculating the variance on each full grid
 *
 */
std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> createBsplineVarianceRefinementOperation(
    size_t degree, size_t numDimensions, sgpp::combigrid::MultiFunction func,
    std::shared_ptr<sgpp::combigrid::LevelManager> levelManager,
    sgpp::combigrid::WeightFunctionsCollection weightFunctions, sgpp::base::DataVector bounds);

/**
 * Creates a level structure according to an averaging level manager using linear calculations on
 * each level as norm. This is a very specific case created for the CO2 example. There it serves as
 * a dummy for storing the levelstructure
 *
 * @param degree        B spline degree
 * @param numDimensions number of dimensions
 * @param func		      the objective function
 * @param levelManager  level manager
 *
 */
std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> createBsplineLinearRefinementOperation(
    size_t degree, size_t numDimensions, sgpp::combigrid::MultiFunction func,
    std::shared_ptr<sgpp::combigrid::LevelManager> levelManager);

/**
 * creates a B spline interpolation operation from a storage of interpolation coefficients
 */
std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> createBsplineLinearCoefficientOperation(
    size_t degree, size_t numDimensions,
    std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> coefficientStorage);

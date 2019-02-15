// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineModifiedGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/datadriven/activeSubspaces/ResponseSurface.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

/**
 * stores a sparse grid not a knot B-spline interpolant in the framework of a repsonse surface
 */
class SparseGridResponseSurfaceNakBspline : public ResponseSurface {
 public:
  /**
   * Constructor
   *
   * @param objectiveFunc		objective Function
   * @param gridType			type of the interpolants grid/basis
   * @param degree			degree of the interpolants basis
   */
  SparseGridResponseSurfaceNakBspline(
      std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc,
      sgpp::base::GridType gridType, size_t degree = 3)
      : ResponseSurface(), objectiveFunc(objectiveFunc), gridType(gridType), degree(degree) {
    initialize();
  }

  /**
   * Destructor
   */
  virtual ~SparseGridResponseSurfaceNakBspline() {}

  /**
   * sets numDim, grid and basis according to objectiveFunction and gridType
   */
  void initialize();

  /**
   * creates a regualar sparse grid interpolant
   * @param level	level of the regular sparse grid
   */
  void createRegularResponseSurface(size_t level);

  /**
   * creates a surplus adaptive sparse grid inteprolant
   * @param maxNumGridPoints	maximum number of grid points of the interpolants grid
   * @param initialLevel		first a regular grid of initialLevel is created.
   * @param refinementsNum		max number of grid points, added in each refinement step
   */
  void createSurplusAdaptiveResponseSurface(size_t maxNumGridPoints, size_t initialLevel,
                                            size_t refinementsNum = 3);

  /**
   * creates a surplus adaptive sparse grid regression approximation
   * @param level				sparse grid level
   * @param evaluationPoints	data points
   * @param functionValues		data values
   * @param lambda				Tikhonov regularization parameter
   */
  void createRegularResponseSurfaceData(size_t level, sgpp::base::DataMatrix evaluationPoints,
                                        sgpp::base::DataVector functionValues,
                                        double lambda = 1e-6);

  /**
   * evaluates this response surface
   * @param v	point to evaluate in
   * @return	the evaluation
   */
  double eval(sgpp::base::DataVector v);

  /**
   * evaluates the response surface and its gradient
   * @param v			point to evaluate in
   * @param gradient	reference to return the repsonse surfaces gradient evaluted in v
   * @return 			the evaluation
   */
  double evalGradient(sgpp::base::DataVector v, sgpp::base::DataVector& gradient);

  /**
   * return the integral of the response surface
   */
  double getIntegral();

  /**
   * @return the interpolation grid
   */
  std::shared_ptr<sgpp::base::Grid> getGrid() { return grid; }

  /**
   * @return the interpolation coefficients
   */
  sgpp::base::DataVector getCoefficients() { return coefficients; }

 private:
  // objective function
  std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc;
  // type of grid/basis
  sgpp::base::GridType gridType;
  // degree of the basis
  size_t degree;
  // number of dimensions
  size_t numDim;
  // the interpolation grid
  std::shared_ptr<sgpp::base::Grid> grid;
  // the interpolation basis
  std::unique_ptr<sgpp::base::SBasis> basis;
  // the interpolation coefficients
  sgpp::base::DataVector coefficients;

  /**
   *refines the grid surplus adaptive and recalculates the interpoaltion coefficients
   *@param refinementsNum	number of grid points which should be refined
   */
  void refineSurplusAdaptive(size_t refinementsNum);

  /**
   * calculates the interpolation coefficients on a given grid
   */
  void calculateInterpolationCoefficients();
};

}  // namespace datadriven
}  // namespace sgpp

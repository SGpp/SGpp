// Copyright(C) 2008 - today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/grid/LevelIndexTypes.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/grid/type/NakPBsplineGrid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/base/tools/Distribution.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/solver/Armadillo.hpp>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/solver/Eigen.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/function/scalar/ResponseSurface.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/optimizer/unconstrained/GradientDescent.hpp>

namespace sgpp {
namespace optimization {

/**
 * stores a sparse grid not a knot B-spline interpolant in the framework of a respsonse surface
 */
class SplineResponseSurface : public ResponseSurface {
 public:
  /**
   * Constructor
   *
   * @param objectiveFunc		objective Function
   * @param gridType			type of the interpolants grid/basis
   * @param degree				degree of the interpolants basis
   */
  SplineResponseSurface(std::shared_ptr<sgpp::base::ScalarFunction> objectiveFunc,
                        sgpp::base::DataVector lb, sgpp::base::DataVector ub,
                        sgpp::base::GridType gridType, size_t degree = 3)
      : ResponseSurface(objectiveFunc->getNumberOfParameters()),
        objectiveFunc(objectiveFunc),
        gridType(gridType),
        degree(degree) {
    this->lb = lb;
    this->ub = ub;
    // dummy values for mean and variance
    mean = 777;
    variance = -1;
    computedMeanFlag = false;
    unitLBounds = sgpp::base::DataVector(numDim, 0.0);
    unitUBounds = sgpp::base::DataVector(numDim, 1.0);
    if (gridType == sgpp::base::GridType::Bspline) {
      grid = std::make_shared<sgpp::base::BsplineGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SBsplineBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::BsplineBoundary) {
      grid = std::make_shared<sgpp::base::BsplineBoundaryGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SBsplineBoundaryBase>(degree);
      boundary = true;
    } else if (gridType == sgpp::base::GridType::ModBspline) {
      grid = std::make_shared<sgpp::base::ModBsplineGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SBsplineModifiedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::BsplineClenshawCurtis) {
      grid = std::make_shared<sgpp::base::BsplineClenshawCurtisGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SBsplineClenshawCurtisBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::FundamentalSpline) {
      grid = std::make_shared<sgpp::base::FundamentalSplineGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SFundamentalSplineBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::ModFundamentalSpline) {
      grid = std::make_shared<sgpp::base::ModFundamentalSplineGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SFundamentalSplineModifiedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakBspline) {
      grid = std::make_shared<sgpp::base::NakBsplineGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SNakBsplineBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
      grid = std::make_shared<sgpp::base::NakBsplineBoundaryGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SNakBsplineBoundaryBase>(degree);
      boundary = true;
    } else if (gridType == sgpp::base::GridType::ModNakBspline) {
      grid = std::make_shared<sgpp::base::ModNakBsplineGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
      grid = std::make_shared<sgpp::base::NakBsplineExtendedGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SNakBsplineExtendedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakPBspline) {
      grid = std::make_shared<sgpp::base::NakPBsplineGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SNakPBsplineBase>(degree);
      boundary = false;
    } else {
      throw sgpp::base::generation_exception("SplineResponseSurface: gridType not supported.");
    }
  }

  /**
   * Constructor
   *
   * @param objectiveFunc		objective Function
   * @param gridType			type of the interpolants grid/basis
   * @param degree				degree of the interpolants basis
   */
  SplineResponseSurface(std::shared_ptr<sgpp::base::ScalarFunction> objectiveFunc,
                        std::shared_ptr<sgpp::base::Grid> grid, sgpp::base::DataVector lb,
                        sgpp::base::DataVector ub, size_t degree = 3)
      : ResponseSurface(objectiveFunc->getNumberOfParameters()),
        objectiveFunc(objectiveFunc),
        degree(degree),
        grid(grid) {
    this->lb = lb;
    this->ub = ub;
    // dummy values for mean and variance
    mean = 777;
    variance = -1;
    computedMeanFlag = false;
    unitLBounds = sgpp::base::DataVector(numDim, 0.0);
    unitUBounds = sgpp::base::DataVector(numDim, 1.0);
    gridType = grid->getType();
    if (gridType == sgpp::base::GridType::Bspline) {
      basis = std::make_unique<sgpp::base::SBsplineBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::BsplineBoundary) {
      basis = std::make_unique<sgpp::base::SBsplineBoundaryBase>(degree);
      boundary = true;
    } else if (gridType == sgpp::base::GridType::ModBspline) {
      basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::BsplineClenshawCurtis) {
      basis = std::make_unique<sgpp::base::SBsplineClenshawCurtisBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::FundamentalSpline) {
      basis = std::make_unique<sgpp::base::SFundamentalSplineBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::ModFundamentalSpline) {
      basis = std::make_unique<sgpp::base::SFundamentalSplineModifiedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakBspline) {
      basis = std::make_unique<sgpp::base::SNakBsplineBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
      basis = std::make_unique<sgpp::base::SNakBsplineBoundaryBase>(degree);
      boundary = true;
    } else if (gridType == sgpp::base::GridType::ModNakBspline) {
      basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
      basis = std::make_unique<sgpp::base::SNakBsplineExtendedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakPBspline) {
      basis = std::make_unique<sgpp::base::SNakPBsplineBase>(degree);
      boundary = false;
    } else {
      throw sgpp::base::generation_exception("SplineResponseSurface: gridType not supported.");
    }
    calculateInterpolationCoefficients();
    interpolant =
        std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
    interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
        *grid, coefficients);
  }

  /**
   * Constructor
   *
   * @param objectiveFunc		objective Function
   * @param gridType			type of the interpolants grid/basis
   * @param degree				degree of the interpolants basis
   */
  SplineResponseSurface(std::shared_ptr<sgpp::base::Grid> grid, sgpp::base::DataVector coefficients,
                        sgpp::base::DataVector lb, sgpp::base::DataVector ub, size_t degree = 3)
      : ResponseSurface(grid->getDimension()),
        degree(degree),
        grid(grid),
        coefficients(coefficients) {
    this->lb = lb;
    this->ub = ub;
    // dummy values for mean and variance
    mean = 777;
    variance = -1;
    computedMeanFlag = false;
    unitLBounds = sgpp::base::DataVector(numDim, 0.0);
    unitUBounds = sgpp::base::DataVector(numDim, 1.0);
    gridType = grid->getType();
    if (gridType == sgpp::base::GridType::Bspline) {
      basis = std::make_unique<sgpp::base::SBsplineBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::BsplineBoundary) {
      basis = std::make_unique<sgpp::base::SBsplineBoundaryBase>(degree);
      boundary = true;
    } else if (gridType == sgpp::base::GridType::ModBspline) {
      basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::BsplineClenshawCurtis) {
      basis = std::make_unique<sgpp::base::SBsplineClenshawCurtisBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::FundamentalSpline) {
      basis = std::make_unique<sgpp::base::SFundamentalSplineBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::ModFundamentalSpline) {
      basis = std::make_unique<sgpp::base::SFundamentalSplineModifiedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakBspline) {
      basis = std::make_unique<sgpp::base::SNakBsplineBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
      basis = std::make_unique<sgpp::base::SNakBsplineBoundaryBase>(degree);
      boundary = true;
    } else if (gridType == sgpp::base::GridType::ModNakBspline) {
      basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
      basis = std::make_unique<sgpp::base::SNakBsplineExtendedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakPBspline) {
      basis = std::make_unique<sgpp::base::SNakPBsplineBase>(degree);
      boundary = false;
    } else {
      throw sgpp::base::generation_exception("SplineResponseSurface: gridType not supported.");
    }
    interpolant =
        std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
    interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
        *grid, coefficients);
  }

  /**
   * Destructor
   */
  ~SplineResponseSurface() {}

  /**
   * creates a regualar sparse grid interpolant
   * @param level	level of the regular sparse grid
   */
  void regular(size_t level);

  /**
   * creates a regular grid of the lowest level with more grid points than specified
   *
   * @param numPoints	desired number of grid points
   * @param verbose		print extra info
   */
  void regularByPoints(size_t numPoints, bool verbose = false);

  /**
   * creates a surplus adaptive sparse grid inteprolant
   * @param maxNumGridPoints	maximum number of grid points of the interpolants grid
   * @param initialLevel		first a regular grid of initialLevel is created.
   * @param refinementsNum		max number of grid points, added in each refinement step
   * @param verbose		print extra info
   */
  void surplusAdaptive(size_t maxNumGridPoints, size_t initialLevel, size_t refinementsNum = 3,
                       bool verbose = false);

  /**
   * creates an adaptive grid based on Ritter-Novak
   * this is favourable for optimization
   *
   * @param maxNumGridPoints	maximum number of grid points of the interpolants grid
   * @param gamma 				    Ritter Novak adaptivity parameter between 0 and1
   * @param initialLevel      initial level for the refinement
   * @param verbose				    print extra info
   */
  void ritterNovak(size_t maxNumGridPoints, double gamma, size_t initialLevel = 3,
                   bool verbose = false);

  /**
   * creates a surplus adaptive sparse grid regression approximation
   * @param level				sparse grid level
   * @param evaluationPoints	data points
   * @param functionValues		data values
   * @param lambda				Tikhonov regularization parameter
   */
  // void regularData(size_t level, sgpp::base::DataMatrix evaluationPoints,
  //                 sgpp::base::DataVector functionValues, double lambda = 1e-6);

  /**
   * evaluates this response surface
   * @param v	point to evaluate in
   * @return	the evaluation
   */
  double eval(sgpp::base::DataVector v) override;

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
   * return the mean of the response surface w.r.t. a probability density function
   *
   * @param 	pdf			the probability density function
   * @param     quadOrder	order of the Gauss Legendre quadrature
   */
  double getMean(sgpp::base::DistributionsVector pdfs, size_t quadOrder);
  /**
   * return the variance of the response surface w.r.t. a probability density function
   *
   * @param 	pdf			the probability density function
   * @param     quadOrder	order of the Gauss Legendre quadrature
   *
   * @return	vector [v,mS,m] containing the variance v and the meanSquare mS and mean m used to
   * compute it
   */
  sgpp::base::DataVector getVariance(sgpp::base::DistributionsVector pdfs, size_t quadOrder);

  /**
   * calculates the minimimum with gradient descent
   *
   * @return argmin of the response surface
   */
  sgpp::base::DataVector optimize();

  /**
   * @return the interpolation grid
   */
  std::shared_ptr<sgpp::base::Grid> getGrid() { return grid; }

  /**
   * @return the interpolation coefficients
   */
  sgpp::base::DataVector getCoefficients() { return coefficients; }

  /**
   * @return the number of grid points
   */
  size_t getSize() { return grid->getSize(); }

  /**
   * @return the lower bounds of the domain
   */
  sgpp::base::DataVector getLowerBounds() { return lb; }

  /**
   * @return the upper bounds of the domain
   */
  sgpp::base::DataVector getUpperBounds() { return ub; }

 private:
  // objective function
  std::shared_ptr<sgpp::base::ScalarFunction> objectiveFunc;
  // type of grid/basis
  sgpp::base::GridType gridType;
  // degree of the basis
  size_t degree;
  // the interpolation grid
  std::shared_ptr<sgpp::base::Grid> grid;
  // the interpolation basis
  std::shared_ptr<sgpp::base::SBasis> basis;
  // the interpolation coefficients
  sgpp::base::DataVector coefficients;
  // function vlaues in the grid points
  sgpp::base::DataVector functionValues;
  // whether or not the grid includes the boundary
  bool boundary;
  // mean value
  double mean;
  // variance value
  double variance;
  // mean computation flag for variance computation
  bool computedMeanFlag;
  sgpp::base::DataVector unitLBounds;
  sgpp::base::DataVector unitUBounds;

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

}  // namespace optimization
}  // namespace sgpp

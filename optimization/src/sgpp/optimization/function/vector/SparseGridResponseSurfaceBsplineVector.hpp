// Copyright(C) 2008 - today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include <sgpp/base/function/vector/InterpolantVectorFunction.hpp>
#include <sgpp/base/function/vector/InterpolantVectorFunctionGradient.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/grid/generation/functors/VectorSurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineModifiedGrid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/base/tools/Distribution.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/solver/Armadillo.hpp>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/function/vector/ResponseSurfaceVector.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/optimizer/unconstrained/GradientDescent.hpp>

namespace sgpp {
namespace optimization {

/**
 * stores a sparse grid not a knot B-spline interpolant in the framework of a respsonse surface
 * for a vector valued objective function (in particular time dependent objective functions, where
 * the return values are the time series)
 * todo (rehmemk): THINK OF A ADAPTIVE REFINEMENT CRITERION! CURRENTLY ONLY REGULAR IMPLEMENTED.
 */
class SparseGridResponseSurfaceBsplineVector : public ResponseSurfaceVector {
 public:
  /**
   * Constructor
   *
   * @param objectiveFunc		objective Function
   * @param gridType			type of the interpolants grid/basis
   * @param degree				degree of the interpolants basis
   */
  SparseGridResponseSurfaceBsplineVector(std::shared_ptr<sgpp::base::VectorFunction> objectiveFunc,
                                         sgpp::base::DataVector lb, sgpp::base::DataVector ub,
                                         sgpp::base::GridType gridType, size_t degree = 3)
      : ResponseSurfaceVector(objectiveFunc->getNumberOfParameters(),
                              objectiveFunc->getNumberOfComponents()),
        objectiveFunc(objectiveFunc),
        gridType(gridType),
        degree(degree) {
    this->lb = lb;
    this->ub = ub;
    // dummy values for mean and variance
    means.resize(numRes, 777);
    variances.resize(numRes, -1);
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
    } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
      grid = std::make_shared<sgpp::base::NakBsplineModifiedGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
      boundary = false;
    } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
      grid = std::make_shared<sgpp::base::NakBsplineExtendedGrid>(numDim, degree);
      basis = std::make_unique<sgpp::base::SNakBsplineExtendedBase>(degree);
      boundary = false;
    } else {
      throw sgpp::base::generation_exception(
          "SparseGridResponseSurfaceBsplineVector: gridType not supported.");
    }
  }

  /**
   * Destructor
   */
  ~SparseGridResponseSurfaceBsplineVector() {}

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
   * @param gamma 				Ritter Novak adaptivity parameter between 0 and 1
   * @param verbose				print extra info
   */
  // void ritterNovak(size_t maxNumGridPoints, double gamma, bool verbose = false);

  /**
   * evaluates this response surface
   * @param v	point to evaluate in
   * @return	the evaluation
   */
  sgpp::base::DataVector eval(sgpp::base::DataVector v) override;

  /**
   * evaluates the response surface and its gradient
   * @param v			point to evaluate in
   * @param jacobian	reference to return the repsonse surfaces jacobian matrix evaluted in v
   * @return 			the evaluation
   */
  sgpp::base::DataVector evalJacobian(sgpp::base::DataVector v, sgpp::base::DataMatrix& jacobian);

  /**
   * return the integrals of the response surface
   */
  sgpp::base::DataVector getIntegrals();

  /**
   * return the mean of the response surface w.r.t. a probability density function
   *
   * @param 	pdf			the probability density function
   * @param     quadOrder	order of the Gauss Legendre quadrature
   *
   * @return vector containing the means
   */
  // sgpp::base::DataVector getMeans(sgpp::base::DistributionsVector pdfs, size_t quadOrder);

  /**
   * return the variance of the response surface w.r.t. a probability density function
   *
   * @param 	pdf			    the probability density function
   * @param   quadOrder	  order of the Gauss Legendre quadrature
   * @param   means       reference to return the means used to calculate the variances
   * @param   meanSquares reference to return the meanSquares used to calculate the variances
   *
   * @return	vector containing the variances
   * compute it
   */
  // sgpp::base::DataVector getVariances(sgpp::base::DistributionsVector pdfs, size_t quadOrder,
  //                                     sgpp::base::DataVector& means,
  //                                     sgpp::base::DataVector& meanSquares);

  /**
   * calculates the minimimum with gradient descent
   *
   * @return argmin of the response surface
   */
  // sgpp::base::DataVector optimize();

  /**
   * @return the interpolation grid
   */
  std::shared_ptr<sgpp::base::Grid> getGrid() { return grid; }

  /**
   * @return the interpolation coefficients
   */
  sgpp::base::DataMatrix getCoefficients() { return coefficients; }

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
  std::shared_ptr<sgpp::base::VectorFunction> objectiveFunc;
  // type of grid/basis
  sgpp::base::GridType gridType;
  // degree of the basis
  size_t degree;
  // the interpolation grid
  std::shared_ptr<sgpp::base::Grid> grid;
  // the interpolation basis
  std::shared_ptr<sgpp::base::SBasis> basis;
  // the interpolation coefficients of shape numGridPoints x numRes
  // Entry [i,j] is the interpolation coefficient for the i'th basis function (grid point) and j'th
  // result entry
  // So every column are basically interpolation coefficients for the univariate function mapping to
  // the j'th result
  sgpp::base::DataMatrix coefficients;
  // function vlaues at the grid points
  sgpp::base::DataMatrix functionValues;
  // whether or not the grid includes the boundary
  bool boundary;
  // mean values
  sgpp::base::DataVector means;
  // variance values
  sgpp::base::DataVector variances;
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

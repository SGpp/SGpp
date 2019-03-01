// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineModifiedGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/datadriven/activeSubspaces/ASResponseSurface.hpp>
#include <sgpp/datadriven/activeSubspaces/EigenFunctionalities.hpp>
#include <sgpp/datadriven/activeSubspaces/MSplineBasis.hpp>
#include <sgpp/datadriven/activeSubspaces/MSplineNakBsplineScalarProducts.hpp>
#include <sgpp/datadriven/activeSubspaces/NakBsplineScalarProducts.hpp>
#include <sgpp/datadriven/activeSubspaces/ResponseSurface.hpp>
#include <sgpp/datadriven/application/RegressionLearner.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/tools/HaltonSequence.hpp>
#include <sgpp/datadriven/tools/SobolSequence.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Reduced response surface on the active subspace of an objective function. The response surface
 * is based on not a knot B-spline basis functions.
 */
class ASResponseSurfaceNakBspline : public ASResponseSurface {
 public:
  /**
   * Constructor
   *
   * @param W1		eigenvectors defining the active subspace
   * @param gridType	type of the basis
   * @param degree	degree of the basis
   */
  ASResponseSurfaceNakBspline(Eigen::MatrixXd W1, sgpp::base::GridType gridType, size_t degree = 3)
      : ASResponseSurface(W1), gridType(gridType), degree(degree) {
    initialize();
  }

  /**
   * sets grid and basis according to gridType and W1
   */
  void initialize();

  /**
   * creates a regular grid of the dimension of the reduced space ( = # columns of W1)
   * then performs regression for the B-spline coefficients with the given evaluationPoints and
   * functionValues.
   *
   * @param evaluationPoints	set of points
   * @param functionValues		the objective function evaluated at evaluationPoints
   * @param level				level of the regular interpolant
   * @param lambda				regularization parameter
   */
  void createRegularReducedSurfaceFromData(sgpp::base::DataMatrix evaluationPoints,
                                           sgpp::base::DataVector functionValues, size_t level,
                                           double lambda = 1e-06);

  /**
   * createRegularReducedSurfaceFromData explicitl<y creates the matrix and uses methods of the
   * eigen library for a least squares approximation of the coefficients. This method accomplishes
   * the same using the SG++ datadriven routines, which apply the CG method without explicitly
   * creating the matrix. ToDo(rehmemk) When both methods are fully implemented compare them and
   * delete the worse one.
   *
   * @param evaluationPoints  set of points in the original space of the objective function
   * @param functionValues	the objective function evaluated at evaluationPoints
   */
  void createRegularReducedSurfaceFromData_DataDriven(sgpp::base::DataMatrix evaluationPoints,
                                                      sgpp::base::DataVector functionValues,
                                                      size_t level, double lambda = 1e-06,
                                                      double exponentBase = 0.25);

  /**
   * creates an adaptive grid of the dimension of the reduced space ( = # columns of W1)
   * then performs regression for the B-spline coefficients with the given evaluationPoints and
   * functionValues.
   *
   * @param evaluationPoints	set of points
   * @param functionValues		the objective function evaluated at evaluationPoints
   * @param level				level of the regular interpolant
   * @param lambda				regularization parameter
   */
  void createAdaptiveReducedSurfaceFromData(size_t maxNumGridPoints,
                                            sgpp::base::DataMatrix evaluationPoints,
                                            sgpp::base::DataVector functionValues,
                                            size_t initialLevel = 1, size_t refinementsNum = 3,
                                            double lambda = 1e-06);

  /**
   * creates a regular grid of the dimension of the reduced space ( = # columns of W1)
   * then interpolates in these points. The function values are calculated by using the
   * Moore-Penrose (pseudo-) inverse of W1, pinvW1. So for a point x the right hand side is
   * f(pinvW1 x)
   *
   * @param level				level of the regular interpolant
   * @param objectiveFunction	the objective function
   */
  void createRegularReducedSurfaceWithPseudoInverse(
      size_t level, std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc);

  /**
   * creates a surplus adaptive grid of the dimension of the reduced space ( = # columns of W1)
   * then interpolates in these points. The function values are calculated by using the
   * Moore-Penrose (pseudo-) inverse of W1, pinvW1. So for a point x the right hand side is
   * f(pinvW1 x)
   *
   * @param maxNumGridPoints	maximum number of grid points
   * @param objectiveFunc		the objective function
   * @param initialLevel		initial regular level of the interpolant. From then on the
   * grid is created adaptively.
   */
  void createAdaptiveReducedSurfaceWithPseudoInverse(
      size_t maxNumGridPoints, std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc,
      size_t initialLevel = 1, size_t refinementsNum = 3);

  /**
   * evalauted the reduced response surface in a point v from the original parameter space. I.e. it
   * evaluates the response surface in W1T*v
   *
   * @param v	point in the original parameter space
   * @return	evaluation at the (transformed) point
   */
  double eval(sgpp::base::DataVector v);

  /**
   * evaluates the reduced response surface's gradient in a point v from the original parameter
   * space. I.e. it evaluates the response surface's gradient in W1T*v
   *
   * @param v			point in the original parameter space
   * qparam gradient 	reference to return the gradient evaluated at the (transformed) point
   * @return			evaluation at the (transformed) point
   */
  double evalGradient(sgpp::base::DataVector v, sgpp::base::DataVector& gradient);

  /**
   * evaluates the 1D interpolant along the 1D active subspace
   */
  double eval1D(double x);

  /**
   * Calculate an approximation for the high dimensional integral based on this one dimensional
   * response surface using Schoenbergs theorem, stating, that the volume orthogonal to the one
   * dimensional space is given by an M_spline defined on the orthogonal transformation of a simplex
   * and thus triangulating the hypercube into simplices of which we can calcualte the volume.
   *  (ON POLYA FREQUENCY FUNCTIONS IV: THE FUNDAMENTAL SPLINE FUNCTIONS AND THEIR LIMITS,
   * Schoenberg, Curry 1966)
   *
   * @argument quadOrder		order of the Gauss quadrature
   *
   * @return					approximation for the high dimensional integral over
   * the unit hypercube
   */
  double getSplineBasedIntegral(size_t quadOrder = 7);

  /**
   * Same as getSplineBasedIntegral, but interpolates the M-spline giving the volume with B-splines
   * for faster evaluation. Speed up with very little loss in accuracy
   * @param approxLevel		level of the B-splines
   * @param approxDegree	degree of the B-splines
   * @return					approximation for the high dimensional integral over
   * the unit hypercube
   */
  double getApproximateSplineBasedIntegral(size_t approxLevel = 7, size_t approxDegree = 3);

  /**
   * Calculates the same as getSplineBasedIntegral but instead of the Schoeneberg theorem uses an
   * interpolant that represents he orthogonal volume. (usually this interpolant is created with
   * getApproximateSplineBasedIntegral)
   *
   * @param volGrid			grid the interpolant for the volume is defined on
   * @param volCoefficients	coefficients of the interpolant for the volume
   * @param volDegree		degree of the interpolant for the volume
   *
   * @ return approximation for the high dimensional integral over the
   * unit hypercube
   */
  double getIntegralFromVolumeInterpolant(std::shared_ptr<sgpp::base::Grid> volGrid,
                                          sgpp::base::DataVector volCoefficients, size_t volDegree);

  /**
   * @return the interpolation grid
   */
  std::shared_ptr<sgpp::base::Grid> getGrid() { return grid; }

  /**
   * @return the interpolation coefficients
   */
  sgpp::base::DataVector getCoefficients() { return coefficients; }

  /**
   * @return the boundary of the active subspace
   */
  inline sgpp::base::DataVector getBounds() {
    sgpp::base::DataVector bounds(2, leftBound1D);
    bounds[1] = rightBound1D;
    return bounds;
  }

 private:
  // grid Type of the response surface
  sgpp::base::GridType gridType;
  // degree of the response surfaces basis
  size_t degree;
  // dimension of the active subspace
  size_t activeDim;
  // the grid of the response surface
  std::shared_ptr<sgpp::base::Grid> grid;
  // the coefficients of the response surface
  sgpp::base::DataVector coefficients;
  // the basis of the response surface
  std::shared_ptr<sgpp::base::SBasis> basis;
  // if the active subspace is [rightBound1D, leftBound1D]. These depend on how W1 rotates the unit
  // hypercube.
  double rightBound1D = 1.0;
  double leftBound1D = 0.0;
  // information from eigen regression
  double mse = 0;
  sgpp::base::DataVector errorPerBasis;

  // ----------------- auxiliary routines -----------
  /**
   * @return factorial of n (n!)
   */
  int factorial(size_t n);

  /**
   * Triangulates the D-dim unit hypercube into D! simplices.
   * Every permutation pi of {0,..,D} defines one simplex in [0,1]^D via
   * 0<= pi(0)<=...<=pi(D)<=1.
   * These are D! dimplices, tringaulating [0,1]^D.
   * This functions returns the volume of these simplices (its the same for all)
   * and W1^T*V, where V is the matrix of the simplex corners.
   *
   * @param projectedCorners	used to return W1^T*V
   *
   * @return volume of the simplices
   */
  double simplexDecomposition(Eigen::MatrixXd& projectedCorners);

  /**
   * Calculate the volume V(y) of the hyperplane {x \in [0,1]^D : W1^T x = y}, i.e. the hyperplane
   * orthogonal to the active subspace through y. This is calculated using the triangulation of the
   * unit cube and the resulting M-splines => V(y) = Vol(simplex) * \sum_{simplices}
   * M-spline(simplex)
   *
   * @param points				evaluation points. Each entry is one y
   * @param simplexVolume		the volume of the simplices (its the sme for all)
   * @param projectedCorners	the simplex Corners V projected to the ative subspace, i.e. W_1^TV
   */
  sgpp::base::DataVector caclculateVolumeSimplexWise(sgpp::base::DataVector points,
                                                     double simplexVolume,
                                                     Eigen::MatrixXd projectedCorners);

  /**
   *refines the grid surplus adaptive and recalculates the interpoaltion coefficients
   *
   *@param refinementsNum	number of grid points which should be refined
   *@pram objectiveFunc		the objective function
   */
  void refineInterpolationSurplusAdaptive(
      size_t refinementsNum, std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc);

  /**
   * calculates the interpolation coefficients for the reduced response surface on the active
   * subspace for a given grid.
   * For x in the original parameter space y = W1T*x is its according point in the active subspace.
   * Accordingly, to interpolate at points y in the active subsapce we want to evaluate
   * f(inv(W1T)*y). W1 is usually not invertible so we use the least squares approximation
   * x = argmin_x ||W_1^Tx - y||_2 and inteprolate in the paris (y,f(x))
   *
   * @param objectiveFunc	the objective function
   */
  void calculateInterpolationCoefficientsWithPseudoInverse(
      std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc);

  /**
   *refines the grid surplus adaptive based on the coefficients and recalculates the regression
   *coefficients
   *
   *@param refinementsNum	number of grid points which should be refined
   * @param evaluationPoints  set of points in the original space of the objective function
   * @param functionValues	the objective function evaluated at evaluationPoints
   * @param lambda				regularization parameter
   */
  void refineRegressionSurplusAdaptive(size_t refinementsNum,
                                       Eigen::MatrixXd transformedPoints_Eigen,
                                       sgpp::base::DataVector functionValues, double lambda);

  /**
   *refines the grid surplus adaptive based on the residual and recalculates the regression
   *coefficients
   *
   *@param refinementsNum	number of grid points which should be refined
   * @param evaluationPoints  set of points in the original space of the objective function
   * @param functionValues	the objective function evaluated at evaluationPoints
   * @param lambda				regularization parameter
   */
  void refineRegressionErrorAdaptive(size_t refinementsNum, Eigen::MatrixXd transformedPoints_Eigen,
                                     sgpp::base::DataVector functionValues, double lambda);

  /**
   * returns the corners of an arbitrary dimensional hypercube
   * @param dimension	dimension of the hypercube
   * @return 			the corners (column wise)
   */
  Eigen::MatrixXd hypercubeVertices(size_t dimension);

  /**
   * determine the 1D active subspaces bounds by transforming all corners c of the unit hypercube to
   * the active subspace by calculating W1T*c. The minimum and maximum of these values are the
   * active subspaces bounds
   */
  void transformationfor1DActiveSubspace();

  /**
   * Transform evaluation points to active subspace, scale if 1D and cast to Eigen vector type
   * @param evaluationPoints
   *
   * @return	transformed evaluationPoints
   */
  Eigen::MatrixXd prepareDataForEigenRegression(sgpp::base::DataMatrix evaluationPoints);
};

}  // namespace datadriven
}  // namespace sgpp

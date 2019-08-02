// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
// #ifdef USE_EIGEN

#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/datadriven/activeSubspaces/ASMatrixBspline.hpp>
#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Used to create, store and use the matrix C for the detection of active subspaces of an
 * analytically given objective function f using a B-spline interpolant. So C_{i,j} = \int \nabla f
 * \nabla f^T dx \approx \int \nabla \hat{f} \nabla \hat{f}^T dx where \hat{f} is a nak B-spline
 * interpolant for f
 */
class ASMatrixBsplineAnalytic : public ASMatrixBspline {
 public:
  /**
   * Constructor
   *
   * @param objectiveFunction objective Function f
   * @param gridType          type of the grid for the interpolant
   * @param degree            degree for the B-spline basis functions
   */
  ASMatrixBsplineAnalytic(std::shared_ptr<sgpp::base::ScalarFunction> objectiveFunc,
                          sgpp::base::GridType gridType, size_t degree)
      : ASMatrixBspline(objectiveFunc->getNumberOfParameters(), degree, gridType),
        objectiveFunc(objectiveFunc) {}

  /**
   * Create a regular interpolant of the objective function f
   *
   * @param level	level of the underlying grid
   */
  void buildRegularInterpolant(size_t level);

  /**
   * Create a spatially adaptive interpolant of the objective function f
   *
   * @param maxNumGridPoints	upper threshold for the number of grid points
   * @param initialLevel		the refinement needs an initial regular grid of initialLevel
   * @param refinementsNum		maximum number of points refined in one step
   */
  void buildAdaptiveInterpolant(size_t maxNumGridPoints, size_t initialLevel = 1,
                                size_t refinementsNum = 3);

  /**
   * calculates the coefficients for the interpolant based on the objective function and grid.
   * Must be called after every change to the grid!
   */
  void calculateCoefficients();

  /**
   * calculates the l2 error of the interpolant
   *
   * @param numMCPoints number of Monte Carlo points
   *
   * @return l2 error
   */
  double l2InterpolationError(size_t numMCPoints = 10000);

  /**
   * calculates the l2 error of the interpolants gradient
   *
   * @param objectiveFuncGradient	the real gradient function
   * @param numMCPoints             number of Monte Carlo Points
   *
   * @return l2 error of gradient
   */
  sgpp::base::DataVector l2InterpolationGradientError(
      std::shared_ptr<sgpp::base::WrapperScalarFunctionGradient> objectiveFuncGradient,
      size_t numMCPoints = 1000);

 private:
  // objective function
  std::shared_ptr<sgpp::base::ScalarFunction> objectiveFunc;
};

}  // namespace datadriven
}  // namespace sgpp

// #endif /* USE_EIGEN */

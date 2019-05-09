// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
// #ifdef USE_EIGEN

#include <sgpp/base/exception/generation_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/optimization/function/scalar/ASInterpolantScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/ASInterpolantScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include <iostream>

namespace sgpp {
namespace optimization {

/**
 * General response surface. Represents an approximation of some function. Usually the approximation
 * is created via interpolation. (But could also be regression for example)
 */
class ResponseSurface {
 public:
  /**
   * Constructor
   */
  explicit ResponseSurface(size_t numDim) : numDim(numDim) {}

  /**
   * Destructor
   */
  virtual ~ResponseSurface() {}

  /**
   * evaluates the response surface
   * [Attention: this is not alway simply "interpolant->eval". In context of active subspaces, for
   * example, the argument v must first be transformed to the active subspace]
   *
   * @param v	point in which the response surface  shall be evaulated
   * @teurn 	evaluation
   */
  virtual double eval(sgpp::base::DataVector v);

  /**
   * evaluates the response surface and its gradient
   *
   * @param v			point in which the response surface shall be evaluated
   * @param gradient 	reference to return the gradient evaluated in v
   * @return 			evaluation
   */
  virtual double evalGradient(sgpp::base::DataVector v, sgpp::base::DataVector& gradient) = 0;

  /**
   * Calculates the l2 error between interpolant and objective function
   *
   * @param objectiveFunction	the objectiveFunction
   * @param numMCPoints			number of Monte Carlo Points
   *
   * @return 					l2 error
   */
  double l2Error(std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc,
                 size_t numMCPoints = 1000);

  /**
   * Calculates the normalized root mean square error between interpolant and objective function
   *
   * @param objectiveFunction	the objectiveFunction
   * @param numMCPoints			number of Monte Carlo Points
   *
   * @return 					vector [NRMSE, l2 error, min, max]
   */
  sgpp::base::DataVector nrmsError(
      std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc, size_t numMCPoints = 1000);

 protected:
  size_t numDim;
  // lower bounds of the objective function's domain
  sgpp::base::DataVector lb;
  // upper bounds of the objective function's domain
  sgpp::base::DataVector ub;
  std::shared_ptr<sgpp::optimization::ScalarFunction> interpolant;
  std::shared_ptr<sgpp::optimization::ScalarFunctionGradient> interpolantGradient;

  /**
   * transforms a point in hyper-rectangle [lBounds,rBounds] to the hyper-rectangle
   * [newlBounds,newuBounds]
   *
   * @param	v			point in [lBounds,uBounds]
   * @param lBounds		lower bounds
   * @param uBounds		upper bounds
   * @param newlBounds  new lower bounds
   * @param newuBounds  new upper bounds
   */
  void transformPoint(sgpp::base::DataVector& v, sgpp::base::DataVector lBounds,
                      sgpp::base::DataVector uBounds, sgpp::base::DataVector newlBounds,
                      sgpp::base::DataVector newuBounds);
  /**
   * calculates the volume of the tensor product domain given by lb and ub
   *
   * @return volume
   */
  double domainVolume();
};

}  // namespace optimization
}  // namespace sgpp

// #endif /* USE_EIGEN */

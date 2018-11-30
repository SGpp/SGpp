// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
// #ifdef USE_EIGEN

#include <sgpp/base/exception/generation_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/optimization/activeSubspaces/EigenFunctionalities.hpp>
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
  ResponseSurface() {}

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
  virtual double eval(sgpp::base::DataVector v) = 0;

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
   * @return the number of grid points
   */
  size_t getSize() { return interpolant->getSize(); }

 protected:
  std::shared_ptr<sgpp::optimization::ASInterpolantScalarFunction> interpolant;
  std::shared_ptr<sgpp::optimization::ASInterpolantScalarFunctionGradient> interpolantGradient;
};

}  // namespace optimization
}  // namespace sgpp

// #endif /* USE_EIGEN */

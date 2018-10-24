// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
//#ifdef USE_EIGEN

#include <sgpp/base/exception/generation_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/optimization/activeSubspaces/EigenFunctionalities.hpp>
#include <sgpp/optimization/function/scalar/ASInterpolantScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/ASInterpolantScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include <iostream>

namespace sgpp {
namespace optimization {

class ResponseSurface {
 public:
  ResponseSurface(){};

  /**
   * Destructor
   */
  virtual ~ResponseSurface() {}

  virtual double eval(sgpp::base::DataVector v) = 0;
  // evaluates gradient AND function
  virtual double evalGradient(sgpp::base::DataVector v, sgpp::base::DataVector& gradient) = 0;

  double l2Error(sgpp::optimization::WrapperScalarFunction objectiveFunc,
                 size_t numMCPoints = 1000);

  size_t getSize() { return interpolant->getSize(); }

 protected:
  std::unique_ptr<sgpp::optimization::ASInterpolantScalarFunction> interpolant;
  std::unique_ptr<sgpp::optimization::ASInterpolantScalarFunctionGradient> interpolantGradient;
};

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */

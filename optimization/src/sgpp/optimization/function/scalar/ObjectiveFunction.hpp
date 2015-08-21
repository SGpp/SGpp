// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_SCALAR_OBJECTIVEFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_SCALAR_OBJECTIVEFUNCTION_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionHessian.hpp>

namespace SGPP {
  namespace optimization {

    typedef ScalarFunction          ObjectiveFunction;
    typedef ScalarFunctionGradient  ObjectiveGradient;
    typedef ScalarFunctionHessian   ObjectiveHessian;

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_SCALAR_OBJECTIVEFUNCTION_HPP */

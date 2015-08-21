// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_SCALAR_INTERPOLANTOBJECTIVEFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_SCALAR_INTERPOLANTOBJECTIVEFUNCTION_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunctionHessian.hpp>

namespace SGPP {
  namespace optimization {

    typedef InterpolantScalarFunction         InterpolantObjectiveFunction;
    typedef InterpolantScalarFunctionGradient InterpolantObjectiveGradient;
    typedef InterpolantScalarFunctionHessian  InterpolantObjectiveHessian;

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_SCALAR_INTERPOLANTOBJECTIVEFUNCTION_HPP */

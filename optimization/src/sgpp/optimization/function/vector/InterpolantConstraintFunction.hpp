// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_VECTOR_INTERPOLANTCONSTRAINTFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_VECTOR_INTERPOLANTCONSTRAINTFUNCTION_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/vector/InterpolantVectorFunction.hpp>
#include <sgpp/optimization/function/vector/InterpolantVectorFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/InterpolantVectorFunctionHessian.hpp>

namespace SGPP {
  namespace optimization {

    typedef InterpolantVectorFunction         InterpolantConstraintFunction;
    typedef InterpolantVectorFunctionGradient InterpolantConstraintGradient;
    typedef InterpolantVectorFunctionHessian  InterpolantConstraintHessian;

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_VECTOR_INTERPOLANTCONSTRAINTFUNCTION_HPP */

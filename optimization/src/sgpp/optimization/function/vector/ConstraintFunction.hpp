// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_VECTOR_CONSTRAINTFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_VECTOR_CONSTRAINTFUNCTION_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionHessian.hpp>

namespace SGPP {
  namespace optimization {

    typedef VectorFunction          ConstraintFunction;
    typedef VectorFunctionGradient  ConstraintGradient;
    typedef VectorFunctionHessian   ConstraintHessian;

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_VECTOR_CONSTRAINTFUNCTION_HPP */

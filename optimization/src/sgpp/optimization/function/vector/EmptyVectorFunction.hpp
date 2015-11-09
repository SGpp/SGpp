// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTION_HPP

#include <cstddef>
#include <sgpp/optimization/function/vector/WrapperVectorFunction.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Empty implementation of VectorFunction.
     * This is intended as a fill-in for ConstrainedOptimizer, if
     * only equality or inequality constraints are supported.
     */
#if defined _WIN32 && !defined _USE_STATICLIB
    extern __declspec(dllimport) WrapperVectorFunction emptyVectorFunction;
#else
    extern WrapperVectorFunction emptyVectorFunction;
#endif

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTION_HPP */

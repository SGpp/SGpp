// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTIONGRADIENT_HPP
#define SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTIONGRADIENT_HPP

#include <cstddef>
#include <sgpp/optimization/function/vector/WrapperVectorFunctionGradient.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Empty implementation of VectorFunctionGradient.
     * This is intended as a fill-in for ConstrainedOptimizer, if
     * only equality or inequality constraints are supported.
     */
    extern WrapperVectorFunctionGradient emptyVectorFunctionGradient;

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORGRADIENT_HPP */

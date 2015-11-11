// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTION_HPP

#include <sgpp/optimization/function/vector/WrapperVectorFunction.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Singleton containing an empty implementation of VectorFunction.
     * This is intended as a fill-in for ConstrainedOptimizer, if
     * only equality or inequality constraints are supported.
     */
    class EmptyVectorFunction {
      public:
        inline static WrapperVectorFunction& getInstance() {
          static WrapperVectorFunction wrapperVectorFunction(
            0, 0, [](const base::DataVector & x,
          base::DataVector & value) {});
          return wrapperVectorFunction;
        }

      private:
        EmptyVectorFunction() {}
        EmptyVectorFunction(const EmptyVectorFunction&) = delete;
        void operator=(const EmptyVectorFunction&) = delete;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTION_HPP */

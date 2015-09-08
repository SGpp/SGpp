// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTION_HPP

#include <cstddef>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Empty implementation of VectorFunction.
     * This is intended as a fill-in for ConstrainedOptimizer, if
     * only equality or inequality constraints are supported.
     */
    class EmptyVectorFunction : public VectorFunction {
      public:
        /**
         * Constructor.
         */
        EmptyVectorFunction() : VectorFunction(0, 0) {
        }

        virtual ~EmptyVectorFunction() {};

        /**
         * Does nothing.
         *
         * @param[in]  x      ignored
         * @param[out] value  ignored
         */
        void eval(const base::DataVector& x, base::DataVector& value) {}

        /**
         * @param[out] clone pointer to cloned object
         */
        virtual void clone(
          std::unique_ptr<VectorFunction>& clone) const {
          clone = std::unique_ptr<VectorFunction>(
                    new EmptyVectorFunction());
        }
    };

    /// instance of EmptyVectorFunction
    extern EmptyVectorFunction emptyVectorFunction;

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTION_HPP */

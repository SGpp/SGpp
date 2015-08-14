// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYCONSTRAINTFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYCONSTRAINTFUNCTION_HPP

#include <cstddef>
#include <sgpp/optimization/function/vector/ConstraintFunction.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Empty implementation of ConstraintFunction.
     * This is intended as a fill-in for ConstrainedOptimizer, if
     * only equality or inequality constraints are supported.
     */
    class EmptyConstraintFunction : public ConstraintFunction {
      public:
        /**
         * Constructor.
         */
        EmptyConstraintFunction() : ConstraintFunction(0, 0) {
        }

        virtual ~EmptyConstraintFunction() {};

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
                    new EmptyConstraintFunction());
        }
    };

    /// instance of EmptyConstraintFunction
    extern EmptyConstraintFunction emptyConstraintFunction;

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYCONSTRAINTFUNCTION_HPP */

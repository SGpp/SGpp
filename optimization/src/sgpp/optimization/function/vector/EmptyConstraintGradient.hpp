// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYCONSTRAINTGRADIENT_HPP
#define SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYCONSTRAINTGRADIENT_HPP

#include <cstddef>
#include <sgpp/optimization/function/vector/ConstraintFunction.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Empty implementation of ConstraintGradient.
     * This is intended as a fill-in for ConstrainedOptimizer, if
     * only equality or inequality constraints are supported.
     */
    class EmptyConstraintGradient : public ConstraintGradient {
      public:
        /**
         * Constructor.
         */
        EmptyConstraintGradient() : ConstraintGradient(0, 0) {
        }

        virtual ~EmptyConstraintGradient() {};

        /**
         * Does nothing.
         *
         * @param[in]  x          ignored
         * @param[out] value      ignored
         * @param[out] gradient   ignored
         */
        void eval(const base::DataVector& x,
                  base::DataVector& value,
                  base::DataMatrix& gradient) {}

        /**
         * @param[out] clone pointer to cloned object
         */
        virtual void clone(
          std::unique_ptr<VectorFunctionGradient>& clone) const {
          clone = std::unique_ptr<VectorFunctionGradient>(
                    new EmptyConstraintGradient());
        }
    };

    /// instance of EmptyConstraintGradient
    extern EmptyConstraintGradient emptyConstraintGradient;

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYCONSTRAINTGRADIENT_HPP */

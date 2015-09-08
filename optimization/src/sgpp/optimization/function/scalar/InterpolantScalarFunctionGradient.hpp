// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_SCALAR_INTERPOLANTSCALARFUNCTIONGRADIENT_HPP
#define SGPP_OPTIMIZATION_FUNCTION_SCALAR_INTERPOLANTSCALARFUNCTIONGRADIENT_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradient.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Sparse grid interpolant gradient of a scalar-valued function.
     *
     * @see InterpolantScalarFunction
     */
    class InterpolantScalarFunctionGradient : public ScalarFunctionGradient {
      public:
        /**
         * Constructor.
         * Do not destruct the grid before the
         * InterpolantScalarFunctionGradient object!
         *
         * @param grid  sparse grid
         * @param alpha coefficient vector
         */
        InterpolantScalarFunctionGradient(base::Grid& grid,
                                          const base::DataVector& alpha) :
          ScalarFunctionGradient(grid.getStorage()->dim()),
          grid(grid),
          opEvalGradient(op_factory::createOperationNaiveEvalGradient(grid)),
          alpha(alpha) {
        }

        /**
         * Evaluation of the function and its gradient.
         *
         * @param      x        evaluation point \f$\vec{x} \in [0, 1]^d\f$
         * @param[out] gradient gradient
         *                      \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
         * @return              \f$f(\vec{x})\f$
         */
        inline float_t eval(const base::DataVector& x,
                            base::DataVector& gradient) {
          // copy x, necessary due to non-existing const correctness
          // in SGPP::base
          base::DataVector y(x);
          return opEvalGradient->evalGradient(alpha, y, gradient);
        }

        /**
         * @param[out] clone pointer to cloned object
         */
        virtual void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const {
          clone = std::unique_ptr<ScalarFunctionGradient>(
                    new InterpolantScalarFunctionGradient(grid, alpha));
        }

        /**
         * @return coefficient vector
         */
        const base::DataVector& getAlpha() const {
          return alpha;
        }

        /**
         * @param alpha coefficient vector
         */
        void setAlpha(const base::DataVector& alpha) {
          this->alpha = alpha;
        }

      protected:
        /// sparse grid
        base::Grid& grid;
        /// pointer to evaluation operation
        std::unique_ptr<base::OperationNaiveEvalGradient> opEvalGradient;
        /// coefficient vector
        base::DataVector alpha;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_SCALAR_INTERPOLANTSCALARFUNCTIONGRADIENT_HPP */

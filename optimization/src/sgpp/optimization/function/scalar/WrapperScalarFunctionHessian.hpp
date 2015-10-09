// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_SCALAR_WRAPPERSCALARFUNCTIONHESSIAN_HPP
#define SGPP_OPTIMIZATION_FUNCTION_SCALAR_WRAPPERSCALARFUNCTIONHESSIAN_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionHessian.hpp>

#include <cstddef>
#include <memory>
#include <functional>

namespace SGPP {
  namespace optimization {

    /**
     * Implementation of ScalarFunctionHessian that
     * wraps a std::function object.
     */
    class WrapperScalarFunctionHessian : public ScalarFunctionHessian {
      public:
        typedef std::function<float_t(const base::DataVector&,
                                      base::DataVector&,
                                      base::DataMatrix&)>
        FunctionHessianEvalType;

        /**
         * Constructor.
         *
         * @param d         dimension of the domain
         * @param fHessian  function gradient to be wrapped
         */
        WrapperScalarFunctionHessian(size_t d,
                                     FunctionHessianEvalType fHessian) :
          ScalarFunctionHessian(d), fHessian(fHessian) {
        }

        /**
         * @param      x        evaluation point \f$\vec{x} \in [0, 1]^d\f$
         * @param[out] gradient gradient
         *                      \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
         * @param[out] hessian  Hessian matrix
         *                      \f$H_f(\vec{x}) \in \mathbb{R}^{d \times d}\f$
         * @return              \f$f(\vec{x})\f$
         */
        inline float_t eval(const base::DataVector& x,
                            base::DataVector& gradient,
                            base::DataMatrix& hessian) {
          return fHessian(x, gradient, hessian);
        }

        /**
         * @param[out] clone pointer to cloned object
         */
        void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const {
          clone = std::unique_ptr<ScalarFunctionHessian>(
                    new WrapperScalarFunctionHessian(d, fHessian));
        }

      protected:
        /// function Hessian to be wrapped
        FunctionHessianEvalType fHessian;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_SCALAR_WRAPPERSCALARFUNCTIONHESSIAN_HPP */

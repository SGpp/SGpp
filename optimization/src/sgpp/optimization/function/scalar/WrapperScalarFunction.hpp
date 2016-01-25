// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_SCALAR_WRAPPERSCALARFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_SCALAR_WRAPPERSCALARFUNCTION_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>

#include <cstddef>
#include <memory>
#include <functional>

namespace SGPP {
  namespace optimization {

    /**
     * Implementation of ScalarFunction that wraps a std::function object.
     */
    class WrapperScalarFunction : public ScalarFunction {
      public:
        typedef std::function<float_t(const base::DataVector&)>
        FunctionEvalType;

        /**
         * Constructor.
         *
         * @param d         dimension of the domain
         * @param f         function to be wrapped
         */
        WrapperScalarFunction(size_t d, FunctionEvalType f) :
          ScalarFunction(d), f(f) {
        }

        /**
         * @param x     evaluation point \f$\vec{x} \in [0, 1]^d\f$
         * @return      \f$f(\vec{x})\f$
         */
        inline float_t eval(const base::DataVector& x) {
          return f(x);
        }

        /**
         * @param[out] clone pointer to cloned object
         */
        void clone(std::unique_ptr<ScalarFunction>& clone) const {
          clone = std::unique_ptr<ScalarFunction>(
                    new WrapperScalarFunction(d, f));
        }

      protected:
        /// function to be wrapped
        FunctionEvalType f;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_SCALAR_WRAPPERSCALARFUNCTION_HPP */

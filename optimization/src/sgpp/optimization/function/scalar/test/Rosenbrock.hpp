// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_ROSENBROCK_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_ROSENBROCK_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/test/TestFunction.hpp>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Rosenbrock test function.
       *
       * Definition:
       * \f$f(\vec{x}) := \sum_{t=1}^{d-1}
       * \left(100 (x_{t+1} - x_t^2)^2 + (1 - x_t)^2\right)\f$,
       * \f$\vec{x} \in [-5, 10]^d\f$,
       * \f$\vec{x}_{\text{opt}} = (1, \dotsc, 1)^{\mathrm{T}}\f$,
       * \f$f_{\text{opt}} = 0\f$
       * (domain scaled to \f$[0, 1]^d\f$)
       */
      class Rosenbrock : public TestFunction {
        public:
          /**
           * Constructor.
           *
           * @param d     dimension of the domain
           */
          Rosenbrock(size_t d) : TestFunction(d) {
          }

          /**
           * Evaluates the test function.
           *
           * @param x     point \f$\vec{x} \in [0, 1]^d\f$
           * @return      \f$f(\vec{x})\f$
           */
          float_t evalUndisplaced(const base::DataVector& x) {
            float_t result = 0.0;

            float_t xt = 15.0 * x[0] - 5.0;

            for (size_t t = 1; t < d; t++) {
              const float_t xtm1 = xt;
              xt = 15.0 * x[t] - 5.0;

              const float_t tmp1 = xt - xtm1 * xtm1;
              const float_t tmp2 = 1.0 - xtm1;
              result += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
            }

            return result;
          }

          /**
           * Returns minimal point and function value of the test function.
           *
           * @param[out] x minimal point
           *               \f$\vec{x}_{\text{opt}} \in [0, 1]^d\f$
           * @return       minimal function value
           *               \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
           */
          float_t getOptimalPointUndisplaced(base::DataVector& x) {
            x.resize(d);
            x.setAll(0.4);
            return 0.0;
          }

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<ScalarFunction>& clone) const {
            clone = std::unique_ptr<ScalarFunction>(new Rosenbrock(*this));
          }
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_ROSENBROCK_HPP */

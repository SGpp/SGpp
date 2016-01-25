// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_SPHERE_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_SPHERE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/test/TestFunction.hpp>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Sphere test function.
       *
       * Definition:
       * \f$f(\vec{x}) := \lVert \vec{x} \rVert_2^2\f$,
       * \f$\vec{x} \in [-1, 9]^d\f$,
       * \f$\vec{x}_{\text{opt}} = \vec{0}\f$,
       * \f$f_{\text{opt}} = 0\f$
       * (domain scaled to \f$[0, 1]^d\f$)
       */
      class Sphere : public TestFunction {
        public:
          /**
           * Constructor.
           *
           * @param d     dimension of the domain
           */
          Sphere(size_t d) : TestFunction(d) {
          }

          /**
           * Evaluates the test function.
           *
           * @param x     point \f$\vec{x} \in [0, 1]^d\f$
           * @return      \f$f(\vec{x})\f$
           */
          float_t evalUndisplaced(const base::DataVector& x) {
            float_t result = 0.0;

            for (size_t t = 0; t < d; t++) {
              const float_t xt = 10.0 * x[t] - 1.0;
              result += xt * xt;
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
            x.setAll(0.1);
            return 0.0;
          }

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<ScalarFunction>& clone) const {
            clone = std::unique_ptr<ScalarFunction>(new Sphere(*this));
          }
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_SPHERE_HPP */

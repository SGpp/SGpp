// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_EGGHOLDER_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_EGGHOLDER_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/test/TestFunction.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Eggholder test function.
       *
       * Definition:
       * \f$f(\vec{x}) :=
       * -(x_2 + 47) \sin\!\left(\sqrt{\left|x_1/2 + x_2 + 47\right|}\right)
       * - x_1 \sin\!\left(\sqrt{\left|x_1 - x_2 - 47\right|}\right)\f$,
       * \f$\vec{x} \in [-512, 512]^2\f$,
       * \f$\vec{x}_{\text{opt}} = (512, 404.2319)^{\mathrm{T}}\f$,
       * \f$f_{\text{opt}} = -959.6407\f$
       * (domain scaled to \f$[0, 1]^2\f$)
       *
       * The displacement is restricted because the minimal point lies on
       * the boundary of \f$[0, 1]^2\f$.
       */
      class Eggholder : public TestFunction {
        public:
          /**
           * Constructor.
           */
          Eggholder() : TestFunction(2) {
          }

          /**
           * Generate normally distributed pseudorandom displacement with
           * default standard deviation and
           * with the restriction of \f$d_1 = 0\f$.
           */
          void generateDisplacement() {
            generateDisplacement(TestFunction::DEFAULT_STANDARD_DEVIATION);
          }

          /**
           * Generate normally distributed pseudorandom displacement
           * with the restriction of \f$d_1 = 0\f$.
           *
           * @param stdDev standard deviation of the displacement coordinates
           */
          void generateDisplacement(float_t stdDev) {
            TestFunction::generateDisplacement(stdDev);
            displacement[0] = 0.0;
          }

          /**
           * Evaluates the test function.
           *
           * @param x     point \f$\vec{x} \in [0, 1]^2\f$
           * @return      \f$f(\vec{x})\f$
           */
          float_t evalUndisplaced(const base::DataVector& x) {
            const float_t x1 = 1024.0 * x[0] - 512.0;
            const float_t x2 = 1024.0 * x[1] - 512.0;

            return -(x2 + 47.0) *
                   std::sin(std::sqrt(std::abs(x1 / 2.0 + x2 + 47.0))) -
                   x1 * std::sin(std::sqrt(std::abs(x1 - (x2 + 47.0))));
          }

          /**
           * Returns minimal point and function value of the test function.
           *
           * @param[out] x minimal point
           *               \f$\vec{x}_{\text{opt}} \in [0, 1]^2\f$
           * @return       minimal function value
           *               \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
           */
          float_t getOptimalPointUndisplaced(base::DataVector& x) {
            x.resize(2);
            x[0] = 1.0;
            x[1] = 0.8947577;
            return evalUndisplaced(x);
          }

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<ScalarFunction>& clone) const {
            clone = std::unique_ptr<ScalarFunction>(new Eggholder(*this));
          }
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_EGGHOLDER_HPP */

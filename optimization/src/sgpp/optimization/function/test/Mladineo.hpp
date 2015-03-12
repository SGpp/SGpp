// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_MLADINEO_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_MLADINEO_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/test/TestFunction.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Mladineo test function.
       *
       * Definition:
       * \f$f(\vec{x}) := 1 + \frac{1}{2} \lVert \vec{x} \rVert_2^2
       *                    - \cos(10 \ln(2x_1)) \cos(10 \ln(3x_2))\f$,
       * \f$\vec{x} \in [0.01, 1]^2\f$,
       * \f$\vec{x}_{\text{opt}} = (0.0115, 0.0144)^{\mathrm{T}}\}\f$,
       * \f$f_{\text{opt}} = 0.000170\f$
       * (domain scaled to \f$[0, 1]^2\f$)
       *
       * The displacement is restricted because the minimal points lie near
       * the origin.
       */
      class Mladineo : public TestFunction {
        public:
          /**
           * Constructor.
           */
          Mladineo() : TestFunction(2) {
          }

          /**
           * Generate normally distributed pseudorandom displacement with
           * default standard deviation and with the restriction of
           * \f$\vec{d} \in [-0.01, 0] \times [-0.01, 0]\f$.
           */
          void generateDisplacement() {
            generateDisplacement(TestFunction::DEFAULT_STANDARD_DEVIATION);
          }

          /**
           * Generate normally distributed pseudorandom displacement
           * with the restriction of
           * \f$\vec{d} \in [-0.01, 0] \times [-0.01, 0]\f$.
           *
           * @param stdDev standard deviation of the displacement coordinates
           */
          void generateDisplacement(float_t stdDev) {
            do {
              TestFunction::generateDisplacement(stdDev);
            } while ((displacement[0] > 0) || (displacement[0] < -0.01) ||
                     (displacement[1] > 0) || (displacement[1] < -0.01));
          }

          /**
           * Evaluates the test function.
           *
           * @param x     point \f$\vec{x} \in [0, 1]^2\f$
           * @return      \f$f(\vec{x})\f$
           */
          float_t evalUndisplaced(const base::DataVector& x) {
            const float_t x1 = 0.99 * x.get(0) + 0.01;
            const float_t x2 = 0.99 * x.get(1) + 0.01;

            return 1.0 + (x1 * x1 + x2 * x2) / 2.0 -
                   cos(10.0 * log(2.0 * x1)) * cos(10.0 * log(3.0 * x2));
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
            x[0] = 0.001542;
            x[1] = 0.004449;
            return evalUndisplaced(x);
          }

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const {
            clone = std::unique_ptr<ObjectiveFunction>(new Mladineo(*this));
          }
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_MLADINEO_HPP */

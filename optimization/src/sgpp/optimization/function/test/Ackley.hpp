// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_ACKLEY_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_ACKLEY_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/test/Test.hpp>

namespace SGPP {
  namespace optimization {
    namespace function {
      namespace test {

        /**
         * Ackley test function.
         *
         * Definition:
         * \f$f(\vec{x}) := 20 + \mathrm{e}
         *                     - 20 \exp\!\left(-\frac{\lVert \vec{x} \rVert_2}{5\sqrt{d}}\right)
         *                     - \exp\!\left(\frac{1}{d} \sum_{t=1}^d \cos(2\pi x_t)\right)\f$,
         * \f$\vec{x} \in [-1, 9]^d\f$,
         * \f$\vec{x}_{\text{opt}} = \vec{0}\f$,
         * \f$f_{\text{opt}} = 0\f$
         * (domain scaled to \f$[0, 1]^d\f$)
         */
        class Ackley : public Test {
          public:
            /**
             * Constructor.
             *
             * @param d     dimension of the domain
             */
            Ackley(size_t d) : Test(d) {
            }

            /**
             * Evaluates the test function.
             *
             * @param x     point \f$\vec{x} \in [0, 1]^d\f$
             * @return      \f$f(\vec{x})\f$
             */
            float_t evalUndisplaced(const std::vector<float_t>& x) {
              float_t result = 0.0;

              float_t arg1 = 0.0;
              float_t arg2 = 0.0;

              for (size_t t = 0; t < d; t++) {
                const float_t xt = 10.0 * x[t] - 1.0;
                arg1 += xt * xt;
                arg2 += cos(2.0 * M_PI * xt);
              }

              result = 20.0 * (1.0 - std::exp(-0.2 * std::sqrt(arg1 / static_cast<float_t>(d))));
              result += M_E - std::exp(arg2 / static_cast<float_t>(d));

              return result;
            }

            /**
             * Returns minimal point and function value of the test function.
             *
             * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^d\f$
             * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
             */
            float_t getOptimalPointUndisplaced(std::vector<float_t>& x) {
              x = std::vector<float_t>(d, 0.1);
              return 0.0;
            }

            /**
             * @return clone of the object
             */
            virtual Objective* clone() const {
              return new Ackley(*this);
            }
        };

      }
    }
  }
}

#endif

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_GRIEWANK_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_GRIEWANK_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/test/Test.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace function {
      namespace test {

        /**
         * Griewank test function.
         *
         * Definition:
         * \f$f(\vec{x}) := 1 + \frac{\lVert \vec{x} \rVert_2^2}{4000}
         *                    - \prod_{t=1}^d \cos\!\left(\frac{x_t}{\sqrt{t}}\right)\f$,
         * \f$\vec{x} \in [-600, 600]^d\f$,
         * \f$\vec{x}_{\text{opt}} = \vec{0}\f$,
         * \f$f_{\text{opt}} = 0\f$
         * (domain scaled to \f$[0, 1]^d\f$)
         */
        class Griewank : public Test {
          public:
            /**
             * Constructor.
             *
             * @param d     dimension of the domain
             */
            Griewank(size_t d) : Test(d) {
            }

            /**
             * Evaluates the test function.
             *
             * @param x     point \f$\vec{x} \in [0, 1]^d\f$
             * @return      \f$f(\vec{x})\f$
             */
            float_t evalUndisplaced(const std::vector<float_t>& x) {
              float_t result = 1.0;
              float_t tmp = 1.0;

              for (size_t t = 0; t < d; t++) {
                const float_t xt = 1200.0 * x[t] - 600.0;
                result += xt * xt / 4000.0;
                tmp *= cos(xt / sqrt(static_cast<float_t>(t + 1)));
              }

              result -= tmp;
              return result;
            }

            /**
             * Returns minimal point and function value of the test function.
             *
             * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^d\f$
             * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
             */
            float_t getOptimalPointUndisplaced(std::vector<float_t>& x) {
              x = std::vector<float_t>(d, 0.5);
              return 0.0;
            }

            /**
             * @return clone of the object
             */
            virtual Objective* clone() const {
              return new Griewank(*this);
            }
        };

      }
    }
  }
}

#endif

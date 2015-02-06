// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_HIMMELBLAU_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_HIMMELBLAU_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/test/Test.hpp>

namespace SGPP {
  namespace optimization {
    namespace function {
      namespace test {

        /**
         * Himmelblau test function.
         *
         * Definition:
         * \f$f(\vec{x}) := (x_1^2 + x_2 - 11)^2 + (x_1 + x_2^2 - 7)^2\f$,
         * \f$\vec{x} \in [-5, 5]^2\f$,
         * \f$\vec{x}_{\text{opt}} \in \{(3, 2)^{\mathrm{T}}, (-2.8051, 3.1313)^{\mathrm{T}},
         *                         (-3.7793, -3.2832)^{\mathrm{T}}, (3.5844, -1.8481)^{\mathrm{T}}\}\f$,
         * \f$f_{\text{opt}} = 0\f$
         * (domain scaled to \f$[0, 1]^2\f$)
         */
        class Himmelblau : public Test {
          public:
            /**
             * Constructor.
             */
            Himmelblau() : Test(2) {
            }

            /**
             * Evaluates the test function.
             *
             * @param x     point \f$\vec{x} \in [0, 1]^2\f$
             * @return      \f$f(\vec{x})\f$
             */
            float_t evalUndisplaced(const std::vector<float_t>& x) {
              const float_t x1 = 10.0 * x[0] - 5.0;
              const float_t x2 = 10.0 * x[1] - 5.0;

              return (x1 * x1 + x2 - 11.0) * (x1 * x1 + x2 - 11.0) + (x1 + x2 * x2 - 7.0) * (x1 + x2 * x2 - 7.0);
            }

            /**
             * Returns minimal point and function value of the test function.
             *
             * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^2\f$
             * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
             */
            float_t getOptimalPointUndisplaced(std::vector<float_t>& x) {
              x.clear();
              x.push_back(0.8);
              x.push_back(0.7);
              return 0.0;
            }

            /**
             * @param[out] clone pointer to cloned object
             */
            virtual void clone(Objective*& clone) const {
              clone = new Himmelblau(*this);
            }
        };

      }
    }
  }
}

#endif

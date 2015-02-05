// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_SCHWEFEL_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_SCHWEFEL_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/test/Test.hpp>

namespace SGPP {
  namespace optimization {
    namespace function {
      namespace test {

        /**
         * Schwefel test function.
         *
         * Definition:
         * \f$f(\vec{x}) := -\sum_{t=1}^d x_t \sin \sqrt{|x_t|}\f$,
         * \f$\vec{x} \in [-500, 500]^d\f$,
         * \f$\vec{x}_{\text{opt}} = (420.9687, \dotsc, 420.9687)^{\mathrm{T}}\f$,
         * \f$f_{\text{opt}} = -418.9829d\f$
         * (domain scaled to \f$[0, 1]^d\f$)
         */
        class Schwefel : public Test {
          public:
            /**
             * Constructor.
             *
             * @param d     dimension of the domain
             */
            Schwefel(size_t d) : Test(d) {
            }

            /**
             * Evaluates the test function.
             *
             * @param x     point \f$\vec{x} \in [0, 1]^d\f$
             * @return      \f$f(\vec{x})\f$
             */
            float_t evalUndisplaced(const std::vector<float_t>& x) {
              float_t result = 0.0;

              for (size_t t = 0; t < d; t++) {
                const float_t xt = 1000.0 * x[t] - 500.0;
                result -= xt * std::sin(std::sqrt(std::abs(xt)));
              }

              return result;
            }

            /**
             * Returns minimal point and function value of the test function.
             *
             * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^d\f$
             * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
             */
            float_t getOptimalPointUndisplaced(std::vector<float_t>& x) {
              x = std::vector<float_t>(d, 0.920968746359982027311844365);
              return evalUndisplaced(x);
            }

            /**
             * @return clone of the object
             */
            virtual Objective* clone() const {
              return new Schwefel(*this);
            }
        };

      }
    }
  }
}

#endif

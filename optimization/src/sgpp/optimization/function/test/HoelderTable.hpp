// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_HOELDERTABLE_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_HOELDERTABLE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/test/Test.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace function {
      namespace test {

        /**
         * Hoelder table test function.
         *
         * Definition:
         * \f$f(\vec{x}) := -\left|\sin x_1 \cos x_2
         *         \exp\!\left(\left|1 - \frac{\lVert \vec{x} \rVert_2}{\pi}\right|\right)\right|\f$,
         * \f$\vec{x} \in [-10, 10]^2\f$,
         * \f$\vec{x}_{\text{opt}} \in \{(8.0550, \pm 9.6646)^{\mathrm{T}},
         *                               (-8.0550, \pm 9.6646)^{\mathrm{T}}\}\f$,
         * \f$f_{\text{opt}} = -19.2085\f$
         * (domain scaled to \f$[0, 1]^2\f$)
         *
         * The displacement is restricted because the minimal points lie near the corners
         * of \f$[0, 1]^2\f$.
         */
        class HoelderTable : public Test {
          public:
            /**
             * Constructor.
             */
            HoelderTable() : Test(2) {
            }

            /**
             * Generate normally distributed pseudorandom displacement with default standard deviation and
             * with the restriction of \f$\vec{d} \in [-0.005, 0.005] \times [-0.01, 0.01]\f$.
             */
            void generateDisplacement() {
              generateDisplacement(Test::DEFAULT_STANDARD_DEVIATION);
            }

            /**
             * Generate normally distributed pseudorandom displacement
             * with the restriction of \f$\vec{d} \in [-0.005, 0.005] \times [-0.01, 0.01]\f$.
             *
             * @param std_dev   standard deviation of the displacement coordinates
             */
            void generateDisplacement(float_t std_dev) {
              do {
                Test::generateDisplacement(std_dev);
              } while ((displacement[0] > 0.005) || (displacement[0] < -0.005) ||
                       (displacement[1] > 0.01) || (displacement[1] < -0.01));
            }

            /**
             * Evaluates the test function.
             *
             * @param x     point \f$\vec{x} \in [0, 1]^2\f$
             * @return      \f$f(\vec{x})\f$
             */
            float_t evalUndisplaced(const std::vector<float_t>& x) {
              const float_t x1 = 20.0 * x[0] - 10.0;
              const float_t x2 = 20.0 * x[1] - 10.0;

              return -std::abs(std::sin(x1) * std::cos(x2) *
                               std::exp(std::abs(1.0 - std::sqrt(x1*x1 + x2*x2) / M_PI)));
            }

            /**
             * Returns minimal point and function value of the test function.
             *
             * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^2\f$
             * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
             */
            float_t getOptimalPointUndisplaced(std::vector<float_t>& x) {
              x.clear();
              x.push_back(0.902751);
              x.push_back(0.9832295);
              return evalUndisplaced(x);
            }

            /**
             * @param[out] clone pointer to cloned object
             */
            virtual void clone(Objective*& clone) const {
              clone = new HoelderTable(*this);
            }
        };

      }
    }
  }
}

#endif

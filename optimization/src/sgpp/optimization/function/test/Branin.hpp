// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_BRANIN_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_BRANIN_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/test/TestFunction.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Branin test function.
       *
       * Definition:
       * \f$f(\vec{x}) :=
       * \left(x_2 - 5.1 x_1^2/(4\pi^2) + 5 x_1/\pi - 6\right)^2
       * + 10 \left(1 - 1/(8\pi)\right) \cos x_1 + 10\f$,
       * \f$\vec{x} \in [-5, 10] \times [0, 15]\f$,
       * \f$\vec{x}_{\text{opt}} \in
       * \{(-\pi, 12,275)^{\mathrm{T}}, (\pi, 2,275)^{\mathrm{T}},
       * (9,42478, 2,475)^{\mathrm{T}}\}\f$,
       * \f$f_{\text{opt}} = 0.397887\f$
       * (domain scaled to \f$[0, 1]^2\f$)
       */
      class Branin : public TestFunction {
        public:
          /**
           * Constructor.
           */
          Branin() : TestFunction(2) {
          }

          /**
           * Evaluates the test function.
           *
           * @param x     point \f$\vec{x} \in [0, 1]^2\f$
           * @return      \f$f(\vec{x})\f$
           */
          float_t evalUndisplaced(const std::vector<float_t>& x) {
            const float_t x1 = 15.0 * x[0] - 5.0;
            const float_t x2 = 15.0 * x[1];
            const float_t tmp = x2 - 5.1 * x1 * x1 / (4.0 * M_PI * M_PI) +
                                5.0 * x1 / M_PI - 6.0;

            return tmp * tmp +
                   10.0 * (1.0 - 1.0 / (8.0 * M_PI)) * std::cos(x1) + 10.0;
          }

          /**
           * Returns minimal point and function value of the test function.
           *
           * @param[out] x minimal point
           *               \f$\vec{x}_{\text{opt}} \in [0, 1]^2\f$
           * @return       minimal function value
           *               \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
           */
          float_t getOptimalPointUndisplaced(std::vector<float_t>& x) {
            x.clear();
            x.push_back(0.5427728435726528825641);
            x.push_back(0.151666666666666666666666667);
            return evalUndisplaced(x);
          }

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const {
            clone = std::unique_ptr<ObjectiveFunction>(new Branin(*this));
          }
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_BRANIN_HPP */

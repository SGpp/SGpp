// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_MICHALEWICZ_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_MICHALEWICZ_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/test/TestFunction.hpp>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Michalewicz test function.
       *
       * Definition:
       * \f$f(\vec{x}) := -\sin x_1 \sin^{20}\!\left(x_1^2/\pi\right) -
       *                   \sin x_2 \sin^{20}\!\left(2x_2^2/\pi\right)\f$,
       * \f$\vec{x} \in [0, 5]^2\f$,
       * \f$\vec{x}_{\text{opt}} = (2.2029055, \pi/2)^{\mathrm{T}}\}\f$,
       * \f$f_{\text{opt}} = -1.8013\f$
       * (domain scaled to \f$[0, 1]^2\f$)
       */
      class Michalewicz : public TestFunction {
        public:
          /**
           * Constructor.
           */
          Michalewicz() : TestFunction(2) {
          }

          /**
           * Evaluates the test function.
           *
           * @param x     point \f$\vec{x} \in [0, 1]^2\f$
           * @return      \f$f(\vec{x})\f$
           */
          float_t evalUndisplaced(const base::DataVector& x) {
            const float_t x1 = 5.0 * x[0];
            const float_t x2 = 5.0 * x[1];

            return -std::sin(x1) * std::pow(std::sin(x1 * x1 / M_PI), 20.0) -
                   std::sin(x2) * std::pow(std::sin(2.0 * x2 * x2 / M_PI),
                                           20.0);
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
            x[0] = 0.440581104135123915;
            x[1] = M_PI / 10.0;
            return evalUndisplaced(x);
          }

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const {
            clone = std::unique_ptr<ObjectiveFunction>(new Michalewicz(*this));
          }
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_MICHALEWICZ_HPP */

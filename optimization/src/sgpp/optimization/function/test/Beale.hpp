// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_BEALE_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_BEALE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/test/TestFunction.hpp>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Beale test function.
       *
       * Definition:
       * \f$f(\vec{x}) := (1.5 - x_1 (1 - x_2))^2 + (2.25 - x_1 (1 - x_2^2))^2
       *                       + (2.625 - x_1 (1 - x_2^3))^2\f$,
       * \f$\vec{x} \in [-5, 5]^2\f$,
       * \f$\vec{x}_{\text{opt}} = (3, 0.5)^{\mathrm{T}}\f$,
       * \f$f_{\text{opt}} = 0\f$
       * (domain scaled to \f$[0, 1]^2\f$)
       */
      class Beale : public TestFunction {
        public:
          /**
           * Constructor.
           */
          Beale() : TestFunction(2) {
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
            const float_t tmp1 = 1.5 - x1 * (1.0 - x2);
            const float_t tmp2 = 2.25 - x1 * (1.0 - x2 * x2);
            const float_t tmp3 = 2.625 - x1 * (1.0 - x2 * x2 * x2);

            return tmp1 * tmp1 + tmp2 * tmp2 + tmp3 * tmp3;
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
            x.push_back(0.8);
            x.push_back(0.55);
            return 0.0;
          }

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<ObjectiveFunction>& clone) const {
            clone = std::unique_ptr<ObjectiveFunction>(new Beale(*this));
          }
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_BEALE_HPP */

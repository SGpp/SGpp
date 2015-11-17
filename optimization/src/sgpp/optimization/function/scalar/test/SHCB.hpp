// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_SHCB_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_SHCB_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/test/TestFunction.hpp>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * SHCB test function.
       *
       * Definition:
       * \f$f(\vec{x}) := x_1^2 \left(4 - 2.1 x_1^2 + x_1^4/3\right) + x_1 x_2
       *                + 4 x_2^2 \left(x_2^2 - 1\right)\f$,
       * \f$\vec{x} \in [-5, 5]^2\f$,
       * \f$\vec{x}_{\text{opt}} \in \{(0.0898, -0.7127)^{\mathrm{T}},
       *                               (-0.0898, 0.7127)^{\mathrm{T}}\}\f$,
       * \f$f_{\text{opt}} = -1.031628\f$
       * (domain scaled to \f$[0, 1]^2\f$)
       */
      class SHCB : public TestFunction {
        public:
          /**
           * Constructor.
           */
          SHCB();

          /**
           * Destructor.
           */
          virtual ~SHCB() override;

          /**
           * Evaluates the test function.
           *
           * @param x     point \f$\vec{x} \in [0, 1]^2\f$
           * @return      \f$f(\vec{x})\f$
           */
          virtual float_t evalUndisplaced(const base::DataVector& x) override;

          /**
           * Returns minimal point and function value of the test function.
           *
           * @param[out] x minimal point
           *               \f$\vec{x}_{\text{opt}} \in [0, 1]^2\f$
           * @return       minimal function value
           *               \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
           */
          virtual float_t getOptimalPointUndisplaced(base::DataVector& x) override;

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<ScalarFunction>& clone) const override;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_SHCB_HPP */

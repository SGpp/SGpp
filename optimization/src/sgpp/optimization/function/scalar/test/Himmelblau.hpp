// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_HIMMELBLAU_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_HIMMELBLAU_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/test/TestFunction.hpp>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Himmelblau test function.
       *
       * Definition:
       * \f$f(\vec{x}) := (x_1^2 + x_2 - 11)^2 + (x_1 + x_2^2 - 7)^2\f$,
       * \f$\vec{x} \in [-5, 5]^2\f$,
       * \f$\vec{x}_{\text{opt}} \in
       * \{(3, 2)^{\mathrm{T}}, (-2.8051, 3.1313)^{\mathrm{T}},
       * (-3.7793, -3.2832)^{\mathrm{T}}, (3.5844, -1.8481)^{\mathrm{T}}\}\f$,
       * \f$f_{\text{opt}} = 0\f$
       * (domain scaled to \f$[0, 1]^2\f$)
       */
      class Himmelblau : public TestFunction {
        public:
          /**
           * Constructor.
           */
          Himmelblau();

          /**
           * Destructor.
           */
          virtual ~Himmelblau() override;

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

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_HIMMELBLAU_HPP */

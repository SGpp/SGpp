// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_EASOM_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_EASOM_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/test/TestFunction.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Easom test function.
       *
       * Definition:
       * \f$f(\vec{x}) := -\cos x_1 \cos x_2
       * \exp(-(x_1 - \pi)^2 - (x_2 - \pi)^2)\f$,
       * \f$\vec{x} \in [-100, 100]^2\f$,
       * \f$\vec{x}_{\text{opt}} = (\pi, \pi)^{\mathrm{T}}\f$,
       * \f$f_{\text{opt}} = -1\f$
       * (domain scaled to \f$[0, 1]^2\f$)
       */
      class Easom : public TestFunction {
        public:
          /**
           * Constructor.
           */
          Easom();

          /**
           * Destructor.
           */
          virtual ~Easom() override;

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

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_EASOM_HPP */

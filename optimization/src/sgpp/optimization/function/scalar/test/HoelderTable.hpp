// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_HOELDERTABLE_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_HOELDERTABLE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/test/TestFunction.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Hoelder table test function.
       *
       * Definition:
       * \f$f(\vec{x}) := -\left|\sin x_1 \cos x_2
       * \exp\!\left(\left|1 -
       * \frac{\lVert \vec{x} \rVert_2}{\pi}\right|\right)\right|\f$,
       * \f$\vec{x} \in [-10, 10]^2\f$,
       * \f$\vec{x}_{\text{opt}} \in \{(8.0550, \pm 9.6646)^{\mathrm{T}},
       *                               (-8.0550, \pm 9.6646)^{\mathrm{T}}\}\f$,
       * \f$f_{\text{opt}} = -19.2085\f$
       * (domain scaled to \f$[0, 1]^2\f$)
       *
       * The displacement is restricted because the minimal points lie near
       * the corners of \f$[0, 1]^2\f$.
       */
      class HoelderTable : public TestFunction {
        public:
          /**
           * Constructor.
           */
          HoelderTable();

          /**
           * Destructor.
           */
          virtual ~HoelderTable() override;

          /**
           * Generate normally distributed pseudorandom displacement with
           * default standard deviation and with the restriction of
           * \f$\vec{d} \in [-0.005, 0.005] \times [-0.01, 0.01]\f$.
           */
          virtual void generateDisplacement() override;

          /**
           * Generate normally distributed pseudorandom displacement
           * with the restriction of
           * \f$\vec{d} \in [-0.005, 0.005] \times [-0.01, 0.01]\f$.
           *
           * @param stdDev standard deviation of the displacement coordinates
           */
          virtual void generateDisplacement(float_t stdDev) override;

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

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_HOELDERTABLE_HPP */

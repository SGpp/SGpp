// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_SCHWEFEL_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_SCHWEFEL_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/test/TestFunction.hpp>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Schwefel test function.
       *
       * Definition:
       * \f$f(\vec{x}) := -\sum_{t=1}^d x_t \sin \sqrt{|x_t|}\f$,
       * \f$\vec{x} \in [-500, 500]^d\f$,
       * \f$\vec{x}_{\text{opt}} =
       * (420.9687, \dotsc, 420.9687)^{\mathrm{T}}\f$,
       * \f$f_{\text{opt}} = -418.9829d\f$
       * (domain scaled to \f$[0, 1]^d\f$)
       */
      class Schwefel : public TestFunction {
        public:
          /**
           * Constructor.
           *
           * @param d     dimension of the domain
           */
          Schwefel(size_t d);

          /**
           * Destructor.
           */
          virtual ~Schwefel() override;

          /**
           * Evaluates the test function.
           *
           * @param x     point \f$\vec{x} \in [0, 1]^d\f$
           * @return      \f$f(\vec{x})\f$
           */
          virtual float_t evalUndisplaced(const base::DataVector& x) override;

          /**
           * Returns minimal point and function value of the test function.
           *
           * @param[out] x minimal point
           *               \f$\vec{x}_{\text{opt}} \in [0, 1]^d\f$
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

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_SCHWEFEL_HPP */

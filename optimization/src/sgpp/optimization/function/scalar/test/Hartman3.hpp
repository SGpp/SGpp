// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_TEST_HARTMAN3_HPP
#define SGPP_OPTIMIZATION_FUNCTION_TEST_HARTMAN3_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/test/TestFunction.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      /**
       * Hartman3 test function.
       *
       * Definition:
       * \f$f(\vec{x}) := -\sum_{i=1}^4 a_i
       * \exp\!\left(-\sum_{t=1}^3 b_{i,t} (x_t - c_{i,t})^2\right)\f$,
       * \f$\vec{a} = \begin{pmatrix}1\\1.2\\3\\3.2\end{pmatrix}\f$,
       * \f$B :=
       *      \begin{pmatrix}
       *          3 & 10 & 30\\
       *          0.1 & 10 & 35\\
       *          3 & 10 & 30\\
       *          0.1 & 10 & 35
       *      \end{pmatrix}\f$,
       * \f$C :=
       *      \begin{pmatrix}
       *          0.3689 & 0.1170 & 0.2673\\
       *          0.4699 & 0.4387 & 0.7470\\
       *          0.1091 & 0.8732 & 0.5547\\
       *          0.0382 & 0.5743 & 0.8828
       *      \end{pmatrix}\f$,
       * \f$\vec{x} \in [0, 1]^3\f$,
       * \f$\vec{x}_{\text{opt}} =
       * (0.114614, 0.555649, 0.852547)^{\mathrm{T}}\f$,
       * \f$f_{\text{opt}} = -3.862785\f$
       */
      class Hartman3 : public TestFunction {
        public:
          /**
           * Constructor.
           */
          Hartman3();

          /**
           * Destructor.
           */
          virtual ~Hartman3() override;

          /**
           * Evaluates the test function.
           *
           * @param x     point \f$\vec{x} \in [0, 1]^3\f$
           * @return      \f$f(\vec{x})\f$
           */
          virtual float_t evalUndisplaced(const base::DataVector& x) override;

          /**
           * Returns minimal point and function value of the test function.
           *
           * @param[out] x minimal point
           *               \f$\vec{x}_{\text{opt}} \in [0, 1]^3\f$
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

#endif /* SGPP_OPTIMIZATION_FUNCTION_TEST_HARTMAN3_HPP */

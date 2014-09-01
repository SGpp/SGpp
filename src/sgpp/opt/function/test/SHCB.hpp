/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_SHCB_HPP
#define SGPP_OPT_FUNCTION_TEST_SHCB_HPP

#include "opt/function/test/Test.hpp"

namespace sg {
  namespace opt {
    namespace function {
      namespace test {

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
        class SHCB : public Test {
          public:
            /**
             * Constructor.
             */
            SHCB() : Test(2) {
            }

            /**
             * Evaluates the test function.
             *
             * @param x     point \f$\vec{x} \in [0, 1]^2\f$
             * @return      \f$f(\vec{x})\f$
             */
            double evalUndisplaced(const std::vector<double>& x) {
              const double x1 = 10.0 * x[0] - 5.0;
              const double x2 = 10.0 * x[1] - 5.0;

              return x1*x1 * (4.0 - 2.1*x1*x1 + x1*x1*x1*x1/3.0) + x1*x2 + 4.0*x2*x2 * (x2*x2 - 1.0);
            }

            /**
             * Returns minimal point and function value of the test function.
             *
             * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^2\f$
             * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
             */
            double getOptimalPointUndisplaced(std::vector<double>& x) {
              x.clear();
              x.push_back(0.50898420131003180624224905);
              x.push_back(0.42873435969792603666027341858);
              return evalUndisplaced(x);
            }

            /**
             * @return clone of the object
             */
            virtual tools::SmartPointer<Objective> clone() {
              return tools::SmartPointer<Objective>(new SHCB(*this));
            }
        };

      }
    }
  }
}

#endif

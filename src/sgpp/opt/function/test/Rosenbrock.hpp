/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_ROSENBROCK_HPP
#define SGPP_OPT_FUNCTION_TEST_ROSENBROCK_HPP

#include "opt/function/test/Test.hpp"

namespace sg {
  namespace opt {
    namespace function {
      namespace test {

        /**
         * Rosenbrock test function.
         *
         * Definition:
         * \f$f(\vec{x}) := \sum_{t=1}^{d-1} \left(100 (x_{t+1} - x_t^2)^2 + (1 - x_t)^2\right)\f$,
         * \f$\vec{x} \in [-5, 10]^d\f$,
         * \f$\vec{x}_{\text{opt}} = (1, \dotsc, 1)^{\mathrm{T}}\f$,
         * \f$f_{\text{opt}} = 0\f$
         * (domain scaled to \f$[0, 1]^d\f$)
         */
        class Rosenbrock : public Test {
          public:
            /**
             * Constructor.
             *
             * @param d     dimension of the domain
             */
            Rosenbrock(size_t d) : Test(d) {
            }

            /**
             * Evaluates the test function.
             *
             * @param x     point \f$\vec{x} \in [0, 1]^d\f$
             * @return      \f$f(\vec{x})\f$
             */
            double evalUndisplaced(const std::vector<double>& x) {
              double result = 0.0;

              double xt = 15.0 * x[0] - 5.0;

              for (size_t t = 1; t < d; t++) {
                const double xtm1 = xt;
                xt = 15.0 * x[t] - 5.0;

                const double tmp1 = xt - xtm1 * xtm1;
                const double tmp2 = 1.0 - xtm1;
                result += 100.0 * tmp1*tmp1 + tmp2*tmp2;
              }

              return result;
            }

            /**
             * Returns minimal point and function value of the test function.
             *
             * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^d\f$
             * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
             */
            double getOptimalPointUndisplaced(std::vector<double>& x) {
              x = std::vector<double>(d, 0.4);
              return 0.0;
            }

            /**
             * @return clone of the object
             */
            virtual tools::SmartPointer<Objective> clone() {
              return tools::SmartPointer<Objective>(new Rosenbrock(*this));
            }
        };

      }
    }
  }
}

#endif

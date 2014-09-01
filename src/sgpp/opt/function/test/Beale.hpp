/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_BEALE_HPP
#define SGPP_OPT_FUNCTION_TEST_BEALE_HPP

#include "opt/function/test/Test.hpp"

namespace sg {
  namespace opt {
    namespace function {
      namespace test {

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
        class Beale : public Test {
          public:
            /**
             * Constructor.
             */
            Beale() : Test(2) {
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
              const double tmp1 = 1.5 - x1 * (1.0 - x2);
              const double tmp2 = 2.25 - x1 * (1.0 - x2*x2);
              const double tmp3 = 2.625 - x1 * (1.0 - x2*x2*x2);

              return tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;
            }

            /**
             * Returns minimal point and function value of the test function.
             *
             * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^2\f$
             * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
             */
            double getOptimalPointUndisplaced(std::vector<double>& x) {
              x.clear();
              x.push_back(0.8);
              x.push_back(0.55);
              return 0.0;
            }

            /**
             * @return clone of the object
             */
            virtual tools::SmartPointer<Objective> clone() {
              return tools::SmartPointer<Objective>(new Beale(*this));
            }
        };

      }
    }
  }
}

#endif

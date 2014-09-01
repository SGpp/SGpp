/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_BRANIN_HPP
#define SGPP_OPT_FUNCTION_TEST_BRANIN_HPP

#include "opt/function/test/Test.hpp"

#include <cmath>

namespace sg {
  namespace opt {
    namespace function {
      namespace test {

        /**
         * Branin test function.
         *
         * Definition:
         * \f$f(\vec{x}) := \left(x_2 - 5.1 x_1^2/(4\pi^2) + 5 x_1/\pi - 6\right)^2
         *                + 10 \left(1 - 1/(8\pi)\right) \cos x_1 + 10\f$,
         * \f$\vec{x} \in [-5, 10] \times [0, 15]\f$,
         * \f$\vec{x}_{\text{opt}} \in
         * \{(-\pi, 12,275)^{\mathrm{T}}, (\pi, 2,275)^{\mathrm{T}}, (9,42478, 2,475)^{\mathrm{T}}\}\f$,
         * \f$f_{\text{opt}} = 0.397887\f$
         * (domain scaled to \f$[0, 1]^2\f$)
         */
        class Branin : public Test {
          public:
            /**
             * Constructor.
             */
            Branin() : Test(2) {
            }

            /**
             * Evaluates the test function.
             *
             * @param x     point \f$\vec{x} \in [0, 1]^2\f$
             * @return      \f$f(\vec{x})\f$
             */
            double evalUndisplaced(const std::vector<double>& x) {
              const double x1 = 15.0 * x[0] - 5.0;
              const double x2 = 15.0 * x[1];
              const double tmp = x2 - 5.1 * x1 * x1 / (4.0 * M_PI * M_PI) + 5.0 * x1 / M_PI - 6.0;

              return tmp * tmp + 10.0 * (1.0 - 1.0 / (8.0 * M_PI)) * std::cos(x1) + 10.0;
            }

            /**
             * Returns minimal point and function value of the test function.
             *
             * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^2\f$
             * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
             */
            double getOptimalPointUndisplaced(std::vector<double>& x) {
              x.clear();
              x.push_back(0.5427728435726528825641);
              x.push_back(0.151666666666666666666666667);
              return evalUndisplaced(x);
            }

            /**
             * @return clone of the object
             */
            virtual tools::SmartPointer<Objective> clone() {
              return tools::SmartPointer<Objective>(new Branin(*this));
            }
        };

      }
    }
  }
}

#endif

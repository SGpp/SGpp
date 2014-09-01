/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_MICHALEWICZ_HPP
#define SGPP_OPT_FUNCTION_TEST_MICHALEWICZ_HPP

#include "opt/function/test/Test.hpp"

namespace sg {
  namespace opt {
    namespace function {
      namespace test {

        /**
         * Michalewicz test function.
         *
         * Definition:
         * \f$f(\vec{x}) := -\sin x_1 \sin^{20}\!\left(x_1^2/\pi\right) -
         *                   \sin x_2 \sin^{20}\!\left(2x_2^2/\pi\right)\f$,
         * \f$\vec{x} \in [0, 5]^2\f$,
         * \f$\vec{x}_{\text{opt}} = (2.2029055, \pi/2)^{\mathrm{T}}\}\f$,
         * \f$f_{\text{opt}} = -1.8013\f$
         * (domain scaled to \f$[0, 1]^2\f$)
         */
        class Michalewicz : public Test {
          public:
            /**
             * Constructor.
             */
            Michalewicz() : Test(2) {
            }

            /**
             * Evaluates the test function.
             *
             * @param x     point \f$\vec{x} \in [0, 1]^2\f$
             * @return      \f$f(\vec{x})\f$
             */
            double evalUndisplaced(const std::vector<double>& x) {
              const double x1 = 5.0 * x[0];
              const double x2 = 5.0 * x[1];

              return -sin(x1) * std::pow(std::sin(x1*x1 / M_PI), 20.0) -
                     sin(x2) * std::pow(std::sin(2.0*x2*x2 / M_PI), 20.0);
            }

            /**
             * Returns minimal point and function value of the test function.
             *
             * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^2\f$
             * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
             */
            double getOptimalPointUndisplaced(std::vector<double>& x) {
              x.clear();
              x.push_back(0.440581104135123915);
              x.push_back(M_PI / 10.0);
              return evalUndisplaced(x);
            }

            /**
             * @return clone of the object
             */
            virtual tools::SmartPointer<Objective> clone() {
              return tools::SmartPointer<Objective>(new Michalewicz(*this));
            }
        };

      }
    }
  }
}

#endif

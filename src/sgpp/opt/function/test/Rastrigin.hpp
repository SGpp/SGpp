/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_RASTRIGIN_HPP
#define SGPP_OPT_FUNCTION_TEST_RASTRIGIN_HPP

#include "opt/function/test/Test.hpp"

#include <cmath>

namespace sg {
  namespace opt {
    namespace function {
      namespace test {

        /**
         * Rastrigin test function.
         *
         * Definition:
         * \f$f(\vec{x}) := 10d + \lVert \vec{x} \rVert_2^2 - 10 \sum_{t=1}^d \cos(2\pi x_t)\f$,
         * \f$\vec{x} \in [-2, 8]^d\f$,
         * \f$\vec{x}_{\text{opt}} = \vec{0}\f$,
         * \f$f_{\text{opt}} = 0\f$
         * (domain scaled to \f$[0, 1]^d\f$)
         */
        class Rastrigin : public Test {
          public:
            /**
             * Constructor.
             *
             * @param d     dimension of the domain
             */
            Rastrigin(size_t d) : Test(d) {
            }

            /**
             * Evaluates the test function.
             *
             * @param x     point \f$\vec{x} \in [0, 1]^d\f$
             * @return      \f$f(\vec{x})\f$
             */
            double evalUndisplaced(const std::vector<double>& x) {
              double result = 10.0 * static_cast<double>(d);

              for (size_t t = 0; t < d; t++) {
                const double xt = 10.0 * x[t] - 2.0;
                result += xt*xt - 10.0 * cos(2*M_PI*xt);
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
              x = std::vector<double>(d, 0.2);
              return 0.0;
            }

            /**
             * @return clone of the object
             */
            virtual tools::SmartPointer<Objective> clone() {
              return tools::SmartPointer<Objective>(new Rastrigin(*this));
            }
        };

      }
    }
  }
}

#endif

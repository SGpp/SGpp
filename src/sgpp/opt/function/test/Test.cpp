/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/function/test/Test.hpp"
#include "opt/tools/RNG.hpp"

#include <iostream>
#include <cstdlib>
#include <cmath>

namespace sg {
  namespace opt {
    namespace function {
      namespace test {

        const double Test::DEFAULT_STANDARD_DEVIATION = 0.01;

        Test::Test(size_t d) :
          Objective(d),
          displacement(std::vector<double>(d, 0.0)) {
        }

        Test::~Test() {
        }

        double Test::eval(const std::vector<double>& x) {
          // displace vector before evaluation
          std::vector<double> x_displaced = x;
          displaceVector(x_displaced);
          return evalUndisplaced(x_displaced);
        }

        double Test::getOptimalPoint(std::vector<double>& x) {
          // reverse displace optimal point
          double fx = getOptimalPointUndisplaced(x);
          reverseDisplaceVector(x);
          return fx;
        }

        void Test::generateDisplacement() {
          generateDisplacement(DEFAULT_STANDARD_DEVIATION);
        }

        void Test::generateDisplacement(double std_dev) {
          this->std_dev = std_dev;

          displacement = std::vector<double>(d, 0.0);

          for (size_t t = 0; t < d; t++) {
            // every component is normally distributed
            displacement[t] = tools::rng.getGaussianRN(std_dev);
          }
        }

        void Test::displaceVector(std::vector<double>& x) const {
          if (x.size() != d) {
            // one could use exceptions for that...
            x.clear();
            std::cerr << "sg::opt::Test::displaceVector: x doesn't match size\n";
            return;
          }

          for (size_t t = 0; t < d; t++) {
            x[t] += displacement[t];
          }
        }

        void Test::reverseDisplaceVector(std::vector<double>& x) const {
          if (x.size() != d) {
            // one could use exceptions for that...
            x.clear();
            std::cerr << "sg::opt::Test::reverseDisplaceVector: x doesn't match size\n";
            return;
          }

          for (size_t t = 0; t < d; t++) {
            x[t] -= displacement[t];
          }
        }

        double Test::getStandardDeviation() const {
          return std_dev;
        }

        void Test::getDisplacement(std::vector<double>& displacement) const {
          displacement = this->displacement;
        }

      }
    }
  }
}

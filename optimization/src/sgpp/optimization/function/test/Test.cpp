// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/test/Test.hpp>
#include <sgpp/optimization/tools/RNG.hpp>

#include <iostream>
#include <cstdlib>
#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace function {
      namespace test {

        const float_t Test::DEFAULT_STANDARD_DEVIATION = 0.01;

        Test::Test(size_t d) :
          Objective(d),
          displacement(std::vector<float_t>(d, 0.0)) {
        }

        Test::~Test() {
        }

        float_t Test::eval(const std::vector<float_t>& x) {
          // displace vector before evaluation
          std::vector<float_t> x_displaced = x;
          displaceVector(x_displaced);
          return evalUndisplaced(x_displaced);
        }

        float_t Test::getOptimalPoint(std::vector<float_t>& x) {
          // reverse displace optimal point
          float_t fx = getOptimalPointUndisplaced(x);
          reverseDisplaceVector(x);
          return fx;
        }

        void Test::generateDisplacement() {
          generateDisplacement(DEFAULT_STANDARD_DEVIATION);
        }

        void Test::generateDisplacement(float_t std_dev) {
          this->std_dev = std_dev;

          displacement = std::vector<float_t>(d, 0.0);

          for (size_t t = 0; t < d; t++) {
            // every component is normally distributed
            displacement[t] = tools::rng.getGaussianRN(std_dev);
          }
        }

        void Test::displaceVector(std::vector<float_t>& x) const {
          if (x.size() != d) {
            // one could use exceptions for that...
            x.clear();
            std::cerr << "SGPP::optimization::Test::displaceVector: x doesn't match size\n";
            return;
          }

          for (size_t t = 0; t < d; t++) {
            x[t] += displacement[t];
          }
        }

        void Test::reverseDisplaceVector(std::vector<float_t>& x) const {
          if (x.size() != d) {
            // one could use exceptions for that...
            x.clear();
            std::cerr << "SGPP::optimization::Test::reverseDisplaceVector: x doesn't match size\n";
            return;
          }

          for (size_t t = 0; t < d; t++) {
            x[t] -= displacement[t];
          }
        }

        float_t Test::getStandardDeviation() const {
          return std_dev;
        }

        void Test::getDisplacement(std::vector<float_t>& displacement) const {
          displacement = this->displacement;
        }

      }
    }
  }
}

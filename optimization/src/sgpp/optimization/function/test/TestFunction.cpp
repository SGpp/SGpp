// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/test/TestFunction.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include <iostream>
#include <cstdlib>
#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      const float_t TestFunction::DEFAULT_STANDARD_DEVIATION = 0.01;

      TestFunction::TestFunction(size_t d) :
        ObjectiveFunction(d),
        stdDev(0.0),
        displacement(std::vector<float_t>(d, 0.0)) {
      }

      TestFunction::~TestFunction() {
      }

      float_t TestFunction::eval(const std::vector<float_t>& x) {
        // displace vector before evaluation
        std::vector<float_t> x_displaced = x;
        displaceVector(x_displaced);
        return evalUndisplaced(x_displaced);
      }

      float_t TestFunction::getOptimalPoint(std::vector<float_t>& x) {
        // reverse displace optimal point
        float_t fx = getOptimalPointUndisplaced(x);
        reverseDisplaceVector(x);
        return fx;
      }

      void TestFunction::generateDisplacement() {
        generateDisplacement(DEFAULT_STANDARD_DEVIATION);
      }

      void TestFunction::generateDisplacement(float_t stdDev) {
        this->stdDev = stdDev;

        displacement = std::vector<float_t>(d, 0.0);

        for (size_t t = 0; t < d; t++) {
          // every component is normally distributed
          displacement[t] = randomNumberGenerator.getGaussianRN(stdDev);
        }
      }

      void TestFunction::displaceVector(std::vector<float_t>& x) const {
        if (x.size() != d) {
          // one could use exceptions for that...
          x.clear();
          std::cerr << "TestFunction::displaceVector: x doesn't match size\n";
          return;
        }

        for (size_t t = 0; t < d; t++) {
          x[t] += displacement[t];
        }
      }

      void TestFunction::reverseDisplaceVector(std::vector<float_t>& x) const {
        if (x.size() != d) {
          // one could use exceptions for that...
          x.clear();
          std::cerr << "TestFunction::reverseDisplaceVector: "
                    << "x doesn't match size\n";
          return;
        }

        for (size_t t = 0; t < d; t++) {
          x[t] -= displacement[t];
        }
      }

      float_t TestFunction::getStandardDeviation() const {
        return stdDev;
      }

      void TestFunction::getDisplacement(
        std::vector<float_t>& displacement) const {
        displacement = this->displacement;
      }

    }
  }
}

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

      TestFunction::TestFunction(size_t d) :
        ObjectiveFunction(d),
        stdDev(0.0),
        displacement(base::DataVector(d, 0.0)) {
      }

      TestFunction::~TestFunction() {
      }

      float_t TestFunction::eval(const base::DataVector& x) {
        // displace vector before evaluation
        base::DataVector x_displaced(x);
        displaceVector(x_displaced);
        return evalUndisplaced(x_displaced);
      }

      float_t TestFunction::getOptimalPoint(base::DataVector& x) {
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

        for (size_t t = 0; t < d; t++) {
          // every component is normally distributed
          displacement[t] = randomNumberGenerator.getGaussianRN(stdDev);
        }
      }

      void TestFunction::displaceVector(base::DataVector& x) const {
        x.resize(d);

        for (size_t t = 0; t < d; t++) {
          x[t] += displacement.get(t);
        }
      }

      void TestFunction::reverseDisplaceVector(base::DataVector& x) const {
        x.resize(d);

        for (size_t t = 0; t < d; t++) {
          x[t] -= displacement.get(t);
        }
      }

      float_t TestFunction::getStandardDeviation() const {
        return stdDev;
      }

      void TestFunction::getDisplacement(base::DataVector& displacement) const {
        displacement.resize(d);
        displacement = this->displacement;
      }

    }
  }
}

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/test/Mladineo.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      Mladineo::Mladineo() : TestFunction(2) {
      }

      Mladineo::~Mladineo() {
      }

      void Mladineo::generateDisplacement() {
        generateDisplacement(TestFunction::DEFAULT_STANDARD_DEVIATION);
      }

      void Mladineo::generateDisplacement(float_t stdDev) {
        do {
          TestFunction::generateDisplacement(stdDev);
        } while ((displacement[0] > 0) || (displacement[0] < -0.01) ||
                 (displacement[1] > 0) || (displacement[1] < -0.01));
      }

      float_t Mladineo::evalUndisplaced(const base::DataVector& x) {
        const float_t x1 = 0.99 * x[0] + 0.01;
        const float_t x2 = 0.99 * x[1] + 0.01;

        return 1.0 + (x1 * x1 + x2 * x2) / 2.0 -
               std::cos(10.0 * std::log(2.0 * x1)) *
               std::cos(10.0 * std::log(3.0 * x2));
      }

      float_t Mladineo::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(2);
        x[0] = 0.001542;
        x[1] = 0.004449;
        return evalUndisplaced(x);
      }

      void Mladineo::clone(std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(new Mladineo(*this));
      }

    }
  }
}

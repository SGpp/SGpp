// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/test/SHCB.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      SHCB::SHCB() : TestFunction(2) {
      }

      SHCB::~SHCB() {
      }

      float_t SHCB::evalUndisplaced(const base::DataVector& x) {
        const float_t x1 = 10.0 * x[0] - 5.0;
        const float_t x2 = 10.0 * x[1] - 5.0;

        return x1 * x1 * (4.0 - 2.1 * x1 * x1 + x1 * x1 * x1 * x1 / 3.0) +
               x1 * x2 + 4.0 * x2 * x2 * (x2 * x2 - 1.0);
      }

      float_t SHCB::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(2);
        x[0] = 0.50898420131003180624224905;
        x[1] = 0.42873435969792603666027341858;
        return evalUndisplaced(x);
      }

      void SHCB::clone(std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(new SHCB(*this));
      }

    }
  }
}

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/test/Hartman6.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      Hartman6::Hartman6() : TestFunction(6) {
      }

      Hartman6::~Hartman6() {
      }

      float_t Hartman6::evalUndisplaced(const base::DataVector& x) {
        return -1.0 * std::exp(-10.0 * (x[0] - 0.1312) * (x[0] - 0.1312) -
                               3.0 * (x[1] - 0.1696) * (x[1] - 0.1696) -
                               17.0 * (x[2] - 0.5569) * (x[2] - 0.5569) -
                               3.5 * (x[3] - 0.0124) * (x[3] - 0.0124) -
                               1.7 * (x[4] - 0.8283) * (x[4] - 0.8283) -
                               8.0 * (x[5] - 0.5886) * (x[5] - 0.5886)) -
               1.2 * std::exp(-0.05 * (x[0] - 0.2329) * (x[0] - 0.2329) -
                              10.0 * (x[1] - 0.4135) * (x[1] - 0.4135) -
                              17.0 * (x[2] - 0.8307) * (x[2] - 0.8307) -
                              0.1 * (x[3] - 0.3736) * (x[3] - 0.3736) -
                              8.0 * (x[4] - 0.1004) * (x[4] - 0.1004) -
                              14.0 * (x[5] - 0.9991) * (x[5] - 0.9991)) -
               3.0 * std::exp(-3.0 * (x[0] - 0.2348) * (x[0] - 0.2348) -
                              3.5 * (x[1] - 0.1451) * (x[1] - 0.1451) -
                              1.7 * (x[2] - 0.3522) * (x[2] - 0.3522) -
                              10.0 * (x[3] - 0.2883) * (x[3] - 0.2883) -
                              17.0 * (x[4] - 0.3047) * (x[4] - 0.3047) -
                              8.0 * (x[5] - 0.6650) * (x[5] - 0.6650)) -
               3.2 * std::exp(-17.0 * (x[0] - 0.4047) * (x[0] - 0.4047) -
                              8.0 * (x[1] - 0.8828) * (x[1] - 0.8828) -
                              0.05 * (x[2] - 0.8732) * (x[2] - 0.8732) -
                              10.0 * (x[3] - 0.5743) * (x[3] - 0.5743) -
                              0.1 * (x[4] - 0.1091) * (x[4] - 0.1091) -
                              14.0 * (x[5] - 0.0381) * (x[5] - 0.0381));
      }

      float_t Hartman6::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(6);
        x[0] = 0.20169;
        x[1] = 0.150011;
        x[2] = 0.476874;
        x[3] = 0.275332;
        x[4] = 0.311652;
        x[5] = 0.6573;
        return evalUndisplaced(x);
      }

      void Hartman6::clone(std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(new Hartman6(*this));
      }

    }
  }
}

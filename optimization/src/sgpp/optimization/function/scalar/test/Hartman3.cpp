// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/test/Hartman3.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      Hartman3::Hartman3() : TestFunction(3) {
      }

      Hartman3::~Hartman3() {
      }

      float_t Hartman3::evalUndisplaced(const base::DataVector& x) {
        return -1.0 * std::exp(-3.0 * (x[0] - 0.3689) * (x[0] - 0.3689) -
                               10.0 * (x[1] - 0.1170) * (x[1] - 0.1170) -
                               30.0 * (x[2] - 0.2673) * (x[2] - 0.2673)) -
               1.2 * std::exp(-0.1 * (x[0] - 0.4699) * (x[0] - 0.4699) -
                              10.0 * (x[1] - 0.4387) * (x[1] - 0.4387) -
                              35.0 * (x[2] - 0.7470) * (x[2] - 0.7470)) -
               3.0 * std::exp(-3.0 * (x[0] - 0.1091) * (x[0] - 0.1091) -
                              10.0 * (x[1] - 0.8732) * (x[1] - 0.8732) -
                              30.0 * (x[2] - 0.5547) * (x[2] - 0.5547)) -
               3.2 * std::exp(-0.1 * (x[0] - 0.0382) * (x[0] - 0.0382) -
                              10.0 * (x[1] - 0.5743) * (x[1] - 0.5743) -
                              35.0 * (x[2] - 0.8828) * (x[2] - 0.8828));
      }

      float_t Hartman3::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(3);
        x[0] = 0.114614;
        x[1] = 0.555649;
        x[2] = 0.852547;
        return evalUndisplaced(x);
      }

      void Hartman3::clone(std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(new Hartman3(*this));
      }

    }
  }
}

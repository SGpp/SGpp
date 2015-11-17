// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/test/Michalewicz.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      Michalewicz::Michalewicz() : TestFunction(2) {
      }

      Michalewicz::~Michalewicz() {
      }

      float_t Michalewicz::evalUndisplaced(const base::DataVector& x) {
        const float_t x1 = 5.0 * x[0];
        const float_t x2 = 5.0 * x[1];

        return -std::sin(x1) * std::pow(std::sin(x1 * x1 / M_PI), 20.0) -
               std::sin(x2) * std::pow(std::sin(2.0 * x2 * x2 / M_PI),
                                       20.0);
      }

      float_t Michalewicz::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(2);
        x[0] = 0.440581104135123915;
        x[1] = M_PI / 10.0;
        return evalUndisplaced(x);
      }

      void Michalewicz::clone(std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(new Michalewicz(*this));
      }

    }
  }
}

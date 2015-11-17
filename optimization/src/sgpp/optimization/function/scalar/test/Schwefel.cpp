// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/test/Schwefel.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      Schwefel::Schwefel(size_t d) : TestFunction(d) {
      }

      Schwefel::~Schwefel() {
      }

      float_t Schwefel::evalUndisplaced(const base::DataVector& x) {
        float_t result = 0.0;

        for (size_t t = 0; t < d; t++) {
          const float_t xt = 1000.0 * x[t] - 500.0;
          result -= xt * std::sin(std::sqrt(std::abs(xt)));
        }

        return result;
      }

      float_t Schwefel::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(d);
        x.setAll(0.920968746359982027311844365);
        return evalUndisplaced(x);
      }

      void Schwefel::clone(std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(new Schwefel(*this));
      }

    }
  }
}

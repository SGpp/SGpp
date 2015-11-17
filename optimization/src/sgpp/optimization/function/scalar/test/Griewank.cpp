// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/test/Griewank.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      Griewank::Griewank(size_t d) : TestFunction(d) {
      }

      Griewank::~Griewank() {
      }

      float_t Griewank::evalUndisplaced(const base::DataVector& x) {
        float_t result = 1.0;
        float_t tmp = 1.0;

        for (size_t t = 0; t < d; t++) {
          const float_t xt = 1200.0 * x[t] - 600.0;
          result += xt * xt / 4000.0;
          tmp *= std::cos(xt / std::sqrt(static_cast<float_t>(t + 1)));
        }

        result -= tmp;
        return result;
      }

      float_t Griewank::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(d);
        x.setAll(0.5);
        return 0.0;
      }

      void Griewank::clone(std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(new Griewank(*this));
      }

    }
  }
}

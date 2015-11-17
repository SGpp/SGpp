// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/test/Rastrigin.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      Rastrigin::Rastrigin(size_t d) : TestFunction(d) {
      }

      Rastrigin::~Rastrigin() {
      }

      float_t Rastrigin::evalUndisplaced(const base::DataVector& x) {
        float_t result = 10.0 * static_cast<float_t>(d);

        for (size_t t = 0; t < d; t++) {
          const float_t xt = 10.0 * x[t] - 2.0;
          result += xt * xt - 10.0 * std::cos(2 * M_PI * xt);
        }

        return result;
      }

      float_t Rastrigin::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(d);
        x.setAll(0.2);
        return 0.0;
      }

      void Rastrigin::clone(std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(new Rastrigin(*this));
      }

    }
  }
}

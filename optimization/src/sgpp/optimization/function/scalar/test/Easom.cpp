// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/test/Easom.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      Easom::Easom() : TestFunction(2) {
      }

      Easom::~Easom() {
      }

      float_t Easom::evalUndisplaced(const base::DataVector& x) {
        const float_t x1 = 200.0 * x[0] - 100.0;
        const float_t x2 = 200.0 * x[1] - 100.0;

        return -std::cos(x1) * std::cos(x2) *
               std::exp(-((x1 - M_PI) * (x1 - M_PI) +
                          (x2 - M_PI) * (x2 - M_PI)));
      }

      float_t Easom::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(2);
        x.setAll(0.51570796326794896619231);
        return -1.0;
      }

      void Easom::clone(std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(new Easom(*this));
      }

    }
  }
}

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/test/Sphere.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_functions {

      Sphere::Sphere(size_t d) : TestFunction(d) {
      }

      Sphere::~Sphere() {
      }

      float_t Sphere::evalUndisplaced(const base::DataVector& x) {
        float_t result = 0.0;

        for (size_t t = 0; t < d; t++) {
          const float_t xt = 10.0 * x[t] - 1.0;
          result += xt * xt;
        }

        return result;
      }

      float_t Sphere::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(d);
        x.setAll(0.1);
        return 0.0;
      }

      void Sphere::clone(std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(new Sphere(*this));
      }

    }
  }
}

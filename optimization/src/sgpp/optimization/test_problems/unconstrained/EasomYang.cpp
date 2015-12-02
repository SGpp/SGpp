// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/EasomYang.hpp>
#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_problems {

      EasomYang::EasomYang(size_t d) :
        UnconstrainedTestProblem(d),
        f(d) {
      }

      EasomYang::~EasomYang() {
      }

      TestScalarFunction& EasomYang::getObjectiveFunction() {
        return f;
      }

      float_t EasomYang::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(d);
        x.setAll(0.75);
        return -1.0;
      }

      EasomYangObjective::EasomYangObjective(size_t d) :
        TestScalarFunction(d) {
      }

      EasomYangObjective::~EasomYangObjective() {
      }

      float_t EasomYangObjective::evalUndisplaced(
        const base::DataVector& x) {
        float_t sum = 0.0;
        float_t product = 1.0;

        for (size_t t = 0; t < d; t++) {
          const float_t xt = 2.0 * M_PI * (2.0 * x[t] - 1.0);
          const float_t diff = xt - M_PI;
          sum += diff * diff;
          product *= -std::cos(xt);
        }

        return -std::exp(-sum) * product;
      }

      void EasomYangObjective::clone(
        std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(
                  new EasomYangObjective(*this));
      }

    }
  }
}

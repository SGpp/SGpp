// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Perm.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_problems {

      Perm::Perm(size_t d) :
        UnconstrainedTestProblem(d),
        f(d) {
      }

      Perm::~Perm() {
      }

      TestScalarFunction& Perm::getObjectiveFunction() {
        return f;
      }

      float_t Perm::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(d);
        const float_t dDbl = static_cast<float_t>(d);

        for (size_t t = 0; t < d; t++) {
          x[t] = 0.5 * (static_cast<float_t>(t + 1) / dDbl + 1.0);
        }

        return 0.0;
      }

      PermObjective::PermObjective(size_t d) :
        TestScalarFunction(d) {
      }

      PermObjective::~PermObjective() {
      }

      float_t PermObjective::evalUndisplaced(
        const base::DataVector& x) {
        float_t result = 0.0;
        const float_t dDbl = static_cast<float_t>(d);

        for (size_t i = 0; i < d; i++) {
          const float_t iDbl = static_cast<float_t>(i + 1);
          float_t innerSum = 0.0;

          for (size_t t = 0; t < d; t++) {
            const float_t xt = dDbl * (2.0 * x[t] - 1.0);
            const float_t tDbl = static_cast<float_t>(t + 1);

            innerSum += (std::pow(tDbl, iDbl) + 1.0) *
                        (std::pow(xt / tDbl, iDbl) - 1.0);
          }

          result += std::pow(innerSum, 2.0);
        }

        return result;
      }

      void PermObjective::clone(
        std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(
                  new PermObjective(*this));
      }

    }
  }
}

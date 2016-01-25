// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/G12.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_problems {

      G12::G12() :
        ConstrainedTestProblem(3),
        f(),
        g(),
        h() {
      }

      G12::~G12() {
      }

      TestScalarFunction& G12::getObjectiveFunction() {
        return f;
      }

      TestVectorFunction& G12::getInequalityConstraintFunction() {
        return g;
      }

      TestVectorFunction& G12::getEqualityConstraintFunction() {
        return h;
      }

      float_t G12::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(3);
        x.setAll(0.5);
        return -1;
      }



      G12Objective::G12Objective() :
        TestScalarFunction(3) {
      }

      G12Objective::~G12Objective() {
      }

      float_t G12Objective::evalUndisplaced(
        const base::DataVector& x) {
        const float_t x1 = 10.0 * x[0];
        const float_t x2 = 10.0 * x[1];
        const float_t x3 = 10.0 * x[2];

        return (std::pow(x1 - 5.0, 2.0) + std::pow(x2 - 5.0, 2.0) +
                std::pow(x3 - 5.0, 2.0)) / 100.0 - 1.0;
      }

      void G12Objective::clone(
        std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(
                  new G12Objective(*this));
      }



      G12InequalityConstraint::G12InequalityConstraint() :
        TestVectorFunction(3, 1) {
      }

      G12InequalityConstraint::~G12InequalityConstraint() {
      }

      void G12InequalityConstraint::evalUndisplaced(
        const base::DataVector& x,
        base::DataVector& value) {
        const float_t x1 = 10.0 * x[0];
        const float_t x2 = 10.0 * x[1];
        const float_t x3 = 10.0 * x[2];
        float_t result = INFINITY;

        for (float_t y1 = 1.0; y1 <= 9.0; y1++) {
          for (float_t y2 = 1.0; y2 <= 9.0; y2++) {
            for (float_t y3 = 1.0; y3 <= 9.0; y3++) {
              const float_t tmp = std::pow(x1 - y1, 2.0) +
                                  std::pow(x2 - y2, 2.0) + std::pow(x3 - y3, 2.0) - 0.0625;

              if (tmp < result) {
                result = tmp;
              }
            }
          }
        }

        value[0] = result;
      }

      void G12InequalityConstraint::clone(
        std::unique_ptr<VectorFunction>& clone) const {
        clone = std::unique_ptr<VectorFunction>(
                  new G12InequalityConstraint(*this));
      }



      G12EqualityConstraint::G12EqualityConstraint() :
        TestVectorFunction(3, 0) {
      }

      G12EqualityConstraint::~G12EqualityConstraint() {
      }

      void G12EqualityConstraint::evalUndisplaced(
        const base::DataVector& x,
        base::DataVector& value) {
      }

      void G12EqualityConstraint::clone(
        std::unique_ptr<VectorFunction>& clone) const {
        clone = std::unique_ptr<VectorFunction>(
                  new G12EqualityConstraint(*this));
      }

    }
  }
}

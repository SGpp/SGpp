// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/G08.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_problems {

      G08::G08() :
        ConstrainedTestProblem(2),
        f(),
        g(),
        h() {
      }

      G08::~G08() {
      }

      TestScalarFunction& G08::getObjectiveFunction() {
        return f;
      }

      TestVectorFunction& G08::getInequalityConstraintFunction() {
        return g;
      }

      TestVectorFunction& G08::getEqualityConstraintFunction() {
        return h;
      }

      float_t G08::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(2);
        x[0] = 0.1227971;
        x[1] = 0.4245373;
        return -0.095825;
      }



      G08Objective::G08Objective() :
        TestScalarFunction(2) {
      }

      G08Objective::~G08Objective() {
      }

      float_t G08Objective::evalUndisplaced(
        const base::DataVector& x) {
        const float_t x1 = 10.0 * x[0];
        const float_t x2 = 10.0 * x[1];

        return -std::pow(std::sin(2.0 * M_PI * x1), 3.0) *
               std::sin(2.0 * M_PI * x2) / (x1 * x1 * x1 * (x1 + x2));
      }

      void G08Objective::clone(
        std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(
                  new G08Objective(*this));
      }



      G08InequalityConstraint::G08InequalityConstraint() :
        TestVectorFunction(2, 2) {
      }

      G08InequalityConstraint::~G08InequalityConstraint() {
      }

      void G08InequalityConstraint::evalUndisplaced(
        const base::DataVector& x,
        base::DataVector& value) {
        const float_t x1 = 10.0 * x[0];
        const float_t x2 = 10.0 * x[1];

        value[0] = x1 * x1 - x2 + 1.0;
        value[1] = 1.0 - x1 + std::pow(x2 - 4.0, 2.0);
      }

      void G08InequalityConstraint::clone(
        std::unique_ptr<VectorFunction>& clone) const {
        clone = std::unique_ptr<VectorFunction>(
                  new G08InequalityConstraint(*this));
      }



      G08EqualityConstraint::G08EqualityConstraint() :
        TestVectorFunction(2, 0) {
      }

      G08EqualityConstraint::~G08EqualityConstraint() {
      }

      void G08EqualityConstraint::evalUndisplaced(
        const base::DataVector& x,
        base::DataVector& value) {
      }

      void G08EqualityConstraint::clone(
        std::unique_ptr<VectorFunction>& clone) const {
        clone = std::unique_ptr<VectorFunction>(
                  new G08EqualityConstraint(*this));
      }

    }
  }
}

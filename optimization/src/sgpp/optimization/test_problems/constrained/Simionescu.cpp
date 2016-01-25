// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/Simionescu.hpp>

#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace test_problems {

      Simionescu::Simionescu() :
        ConstrainedTestProblem(2),
        f(),
        g(),
        h() {
      }

      Simionescu::~Simionescu() {
      }

      TestScalarFunction& Simionescu::getObjectiveFunction() {
        return f;
      }

      TestVectorFunction& Simionescu::getInequalityConstraintFunction() {
        return g;
      }

      TestVectorFunction& Simionescu::getEqualityConstraintFunction() {
        return h;
      }

      float_t Simionescu::getOptimalPointUndisplaced(base::DataVector& x) {
        x.resize(2);
        x[0] = 0.1605887450304572;
        x[1] = 0.8394112549695428;
        return -0.072;
      }



      SimionescuObjective::SimionescuObjective() :
        TestScalarFunction(2) {
      }

      SimionescuObjective::~SimionescuObjective() {
      }

      float_t SimionescuObjective::evalUndisplaced(
        const base::DataVector& x) {
        const float_t x1 = 2.5 * x[0] - 1.25;
        const float_t x2 = 2.5 * x[1] - 1.25;
        return 0.1 * x1 * x2;
      }

      void SimionescuObjective::clone(
        std::unique_ptr<ScalarFunction>& clone) const {
        clone = std::unique_ptr<ScalarFunction>(
                  new SimionescuObjective(*this));
      }



      SimionescuInequalityConstraint::SimionescuInequalityConstraint() :
        TestVectorFunction(2, 1) {
      }

      SimionescuInequalityConstraint::~SimionescuInequalityConstraint() {
      }

      void SimionescuInequalityConstraint::evalUndisplaced(
        const base::DataVector& x,
        base::DataVector& value) {
        const float_t rT = 1;
        const float_t rS = 0.2;
        const float_t n = 8.0;
        const float_t x1 = 2.5 * x[0] - 1.25;
        const float_t x2 = 2.5 * x[1] - 1.25;
        const float_t tmp = rT + rS * std::cos(n * std::atan(x1 / x2));

        value[0] = x1 * x1 + x2 * x2 - tmp * tmp;
      }

      void SimionescuInequalityConstraint::clone(
        std::unique_ptr<VectorFunction>& clone) const {
        clone = std::unique_ptr<VectorFunction>(
                  new SimionescuInequalityConstraint(*this));
      }



      SimionescuEqualityConstraint::SimionescuEqualityConstraint() :
        TestVectorFunction(2, 0) {
      }

      SimionescuEqualityConstraint::~SimionescuEqualityConstraint() {
      }

      void SimionescuEqualityConstraint::evalUndisplaced(
        const base::DataVector& x,
        base::DataVector& value) {
      }

      void SimionescuEqualityConstraint::clone(
        std::unique_ptr<VectorFunction>& clone) const {
        clone = std::unique_ptr<VectorFunction>(
                  new SimionescuEqualityConstraint(*this));
      }

    }
  }
}

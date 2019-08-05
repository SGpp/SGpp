// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/G08.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

G08::G08() : ConstrainedTestProblem(2), f(), g(), h() {}

G08::~G08() {}

TestScalarFunction& G08::getObjectiveFunction() { return f; }

TestVectorFunction& G08::getInequalityConstraintFunction() { return g; }

TestVectorFunction& G08::getEqualityConstraintFunction() { return h; }

double G08::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.363985679168500;
  x[1] = 0.415124455491333;
  return -0.09582504141804;
}

G08Objective::G08Objective() : TestScalarFunction(2) {}

G08Objective::~G08Objective() {}

double G08Objective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 2.0 * x[0] + 0.5;
  const double x2 = 3.0 * x[1] + 3.0;

  return -std::pow(std::sin(2.0 * M_PI * x1), 3.0) * std::sin(2.0 * M_PI * x2) /
         (x1 * x1 * x1 * (x1 + x2));
}

void G08Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new G08Objective(*this));
}

G08InequalityConstraint::G08InequalityConstraint() : TestVectorFunction(2, 2) {}

G08InequalityConstraint::~G08InequalityConstraint() {}

void G08InequalityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {
  const double x1 = 2.0 * x[0] + 0.5;
  const double x2 = 3.0 * x[1] + 3.0;

  value[0] = x1 * x1 - x2 + 1.0;
  value[1] = 1.0 - x1 + std::pow(x2 - 4.0, 2.0);
}

void G08InequalityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G08InequalityConstraint(*this));
}

G08EqualityConstraint::G08EqualityConstraint() : TestVectorFunction(2, 0) {}

G08EqualityConstraint::~G08EqualityConstraint() {}

void G08EqualityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {}

void G08EqualityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G08EqualityConstraint(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

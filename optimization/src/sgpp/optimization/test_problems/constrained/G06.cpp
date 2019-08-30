// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/G06.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

G06::G06() : ConstrainedTestProblem(2), f(), g(), h() {}

G06::~G06() {}

TestScalarFunction& G06::getObjectiveFunction() { return f; }

TestVectorFunction& G06::getInequalityConstraintFunction() { return g; }

TestVectorFunction& G06::getEqualityConstraintFunction() { return h; }

double G06::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.0125862068965517;
  x[1] = 0.00842960789215478;
  return -6961.81387558014;
}

G06Objective::G06Objective() : TestScalarFunction(2) {}

G06Objective::~G06Objective() {}

double G06Objective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 87.0 * x[0] + 13.0;
  const double x2 = 100.0 * x[1];

  return std::pow(x1 - 10.0, 3.0) + std::pow(x2 - 20.0, 3.0);
}

void G06Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new G06Objective(*this));
}

G06InequalityConstraint::G06InequalityConstraint() : TestVectorFunction(2, 2) {}

G06InequalityConstraint::~G06InequalityConstraint() {}

void G06InequalityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {
  const double x1 = 87.0 * x[0] + 13.0;
  const double x2 = 100.0 * x[1];

  value[0] = -std::pow(x1 - 5.0, 2.0) - std::pow(x2 - 5.0, 2.0) + 100;
  value[1] = std::pow(x1 - 6.0, 2.0) + std::pow(x2 - 5.0, 2.0) - 82.81;
}

void G06InequalityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G06InequalityConstraint(*this));
}

G06EqualityConstraint::G06EqualityConstraint() : TestVectorFunction(2, 0) {}

G06EqualityConstraint::~G06EqualityConstraint() {}

void G06EqualityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {}

void G06EqualityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G06EqualityConstraint(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

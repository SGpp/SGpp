// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/G13.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

G13::G13() : ConstrainedTestProblem(5), f(), g(), h() {}

G13::~G13() {}

TestScalarFunction& G13::getObjectiveFunction() { return f; }

TestVectorFunction& G13::getInequalityConstraintFunction() { return g; }

TestVectorFunction& G13::getEqualityConstraintFunction() { return h; }

double G13::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(5);
  x[0] = 0.126708043478261;
  x[1] = 0.846893260869565;
  x[2] = 0.785507343750000;
  x[3] = 0.380681046875000;
  x[4] = 0.380680468750000;
  return 0.0539498;
}

G13Objective::G13Objective() : TestScalarFunction(5) {}

G13Objective::~G13Objective() {}

double G13Objective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 4.6 * x[0] - 2.3;
  const double x2 = 4.6 * x[1] - 2.3;
  const double x3 = 6.4 * x[2] - 3.2;
  const double x4 = 6.4 * x[3] - 3.2;
  const double x5 = 6.4 * x[4] - 3.2;

  return std::exp(x1 * x2 * x3 * x4 * x5);
}

void G13Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new G13Objective(*this));
}

G13InequalityConstraint::G13InequalityConstraint() : TestVectorFunction(5, 0) {}

G13InequalityConstraint::~G13InequalityConstraint() {}

void G13InequalityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {}

void G13InequalityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G13InequalityConstraint(*this));
}

G13EqualityConstraint::G13EqualityConstraint() : TestVectorFunction(5, 3) {}

G13EqualityConstraint::~G13EqualityConstraint() {}

void G13EqualityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {
  const double x1 = 4.6 * x[0] - 2.3;
  const double x2 = 4.6 * x[1] - 2.3;
  const double x3 = 6.4 * x[2] - 3.2;
  const double x4 = 6.4 * x[3] - 3.2;
  const double x5 = 6.4 * x[4] - 3.2;

  value[0] = x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4 + x5 * x5 - 10.0;
  value[1] = x2 * x3 - 5.0 * x4 * x5;
  value[2] = x1 * x1 * x1 + x2 * x2 * x2 + 1.0;
}

void G13EqualityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G13EqualityConstraint(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

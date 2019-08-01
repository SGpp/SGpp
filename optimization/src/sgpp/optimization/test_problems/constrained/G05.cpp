// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/G05.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

G05::G05() : ConstrainedTestProblem(4), f(), g(), h() {}

G05::~G05() {}

TestScalarFunction& G05::getObjectiveFunction() { return f; }

TestVectorFunction& G05::getInequalityConstraintFunction() { return g; }

TestVectorFunction& G05::getEqualityConstraintFunction() { return h; }

double G05::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(4);
  x[0] = 0.566621083333333;
  x[1] = 0.855055833333333;
  x[2] = 0.608069454545455;
  x[3] = 0.139787636363636;
  return 5126.497478059326;
}

G05Objective::G05Objective() : TestScalarFunction(4) {}

G05Objective::~G05Objective() {}

double G05Objective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 1200.0 * x[0];
  const double x2 = 1200.0 * x[1];
  // const double x3 = 1.1 * x[2] - 0.55;
  // const double x4 = 1.1 * x[3] - 0.55;

  return 3.0 * x1 + 0.000001 * std::pow(x1, 3.0) + 2.0 * x2 + (0.000002 / 3.0) * std::pow(x2, 3.0);
}

void G05Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new G05Objective(*this));
}

G05InequalityConstraint::G05InequalityConstraint() : TestVectorFunction(4, 2) {}

G05InequalityConstraint::~G05InequalityConstraint() {}

void G05InequalityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {
  // const double x1 = 1200.0 * x[0];
  // const double x2 = 1200.0 * x[1];
  const double x3 = 1.1 * x[2] - 0.55;
  const double x4 = 1.1 * x[3] - 0.55;

  value[0] = -x4 + x3 - 0.55;
  value[1] = -x3 + x4 - 0.55;
}

void G05InequalityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G05InequalityConstraint(*this));
}

G05EqualityConstraint::G05EqualityConstraint() : TestVectorFunction(4, 3) {}

G05EqualityConstraint::~G05EqualityConstraint() {}

void G05EqualityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {
  const double x1 = 1200.0 * x[0];
  const double x2 = 1200.0 * x[1];
  const double x3 = 1.1 * x[2] - 0.55;
  const double x4 = 1.1 * x[3] - 0.55;

  value[0] = 1000.0 * std::sin(-x3 - 0.25) + 1000.0 * std::sin(-x4 - 0.25) + 894.8 - x1;
  value[1] = 1000.0 * std::sin(x3 - 0.25) + 1000.0 * std::sin(x3 - x4 - 0.25) + 894.8 - x2;
  value[2] = 1000.0 * std::sin(x4 - 0.25) + 1000.0 * std::sin(x4 - x3 - 0.25) + 1294.8;
}

void G05EqualityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G05EqualityConstraint(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

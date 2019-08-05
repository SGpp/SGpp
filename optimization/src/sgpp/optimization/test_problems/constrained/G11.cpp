// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/G11.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

G11::G11() : ConstrainedTestProblem(2), f(), g(), h() {}

G11::~G11() {}

TestScalarFunction& G11::getObjectiveFunction() { return f; }

TestVectorFunction& G11::getInequalityConstraintFunction() { return g; }

TestVectorFunction& G11::getEqualityConstraintFunction() { return h; }

double G11::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.853553390593274;
  x[1] = 0.75;
  return 0.75;
}

G11Objective::G11Objective() : TestScalarFunction(2) {}

G11Objective::~G11Objective() {}

double G11Objective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 2.0 * x[0] - 1.0;
  const double x2 = 2.0 * x[1] - 1.0;

  return x1 * x1 + std::pow(x2 - 1.0, 2.0);
}

void G11Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new G11Objective(*this));
}

G11InequalityConstraint::G11InequalityConstraint() : TestVectorFunction(2, 0) {}

G11InequalityConstraint::~G11InequalityConstraint() {}

void G11InequalityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {}

void G11InequalityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G11InequalityConstraint(*this));
}

G11EqualityConstraint::G11EqualityConstraint() : TestVectorFunction(2, 1) {}

G11EqualityConstraint::~G11EqualityConstraint() {}

void G11EqualityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {
  const double x1 = 2.0 * x[0] - 1.0;
  const double x2 = 2.0 * x[1] - 1.0;

  value[0] = x2 - x1 * x1;
}

void G11EqualityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G11EqualityConstraint(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

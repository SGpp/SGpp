// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/G03.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

G03::G03(size_t d) : ConstrainedTestProblem(d), f(d), g(d), h(d) {}

G03::~G03() {}

TestScalarFunction& G03::getObjectiveFunction() { return f; }

TestVectorFunction& G03::getInequalityConstraintFunction() { return g; }

TestVectorFunction& G03::getEqualityConstraintFunction() { return h; }

double G03::getOptimalPointUndisplaced(base::DataVector& x) {
  const double dDbl = static_cast<double>(d);
  x.resize(d);
  x.setAll(std::pow(dDbl, -0.5));
  return -std::pow(dDbl, -0.5 * dDbl);
}

G03Objective::G03Objective(size_t d) : TestScalarFunction(d) {}

G03Objective::~G03Objective() {}

double G03Objective::evalUndisplaced(const base::DataVector& x) {
  double result = -1.0;

  for (size_t t = 0; t < d; t++) {
    result *= x[t];
  }

  return result;
}

void G03Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new G03Objective(*this));
}

G03InequalityConstraint::G03InequalityConstraint(size_t d) : TestVectorFunction(d, 0) {}

G03InequalityConstraint::~G03InequalityConstraint() {}

void G03InequalityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {}

void G03InequalityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G03InequalityConstraint(*this));
}

G03EqualityConstraint::G03EqualityConstraint(size_t d) : TestVectorFunction(d, 1) {}

G03EqualityConstraint::~G03EqualityConstraint() {}

void G03EqualityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {
  double result = -1.0;

  for (size_t t = 0; t < d; t++) {
    result += x[t] * x[t];
  }

  value[0] = result;
}

void G03EqualityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G03EqualityConstraint(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

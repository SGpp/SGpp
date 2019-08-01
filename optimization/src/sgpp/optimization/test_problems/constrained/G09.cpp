// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/G09.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

G09::G09() : ConstrainedTestProblem(7), f(), g(), h() {}

G09::~G09() {}

TestScalarFunction& G09::getObjectiveFunction() { return f; }

TestVectorFunction& G09::getInequalityConstraintFunction() { return g; }

TestVectorFunction& G09::getEqualityConstraintFunction() { return h; }

double G09::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(7);
  x[0] = 0.61652495;
  x[1] = 0.59756860;
  x[2] = 0.47612293;
  x[3] = 0.71828630;
  x[4] = 0.46877565;
  x[5] = 0.55190655;
  x[6] = 0.57971135;
  return 680.630111240756;
}

G09Objective::G09Objective() : TestScalarFunction(7) {}

G09Objective::~G09Objective() {}

double G09Objective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 20.0 * x[0] - 10.0;
  const double x2 = 20.0 * x[1] - 10.0;
  const double x3 = 20.0 * x[2] - 10.0;
  const double x4 = 20.0 * x[3] - 10.0;
  const double x5 = 20.0 * x[4] - 10.0;
  const double x6 = 20.0 * x[5] - 10.0;
  const double x7 = 20.0 * x[6] - 10.0;

  return std::pow(x1 - 10.0, 2.0) + 5.0 * std::pow(x2 - 12.0, 2.0) + std::pow(x3, 4.0) +
         3.0 * std::pow(x4 - 11.0, 2.0) + 10.0 * std::pow(x5, 6.0) + 7.0 * x6 * x6 +
         std::pow(x7, 4.0) - 4.0 * x6 * x7 - 10.0 * x6 - 8.0 * x7;
}

void G09Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new G09Objective(*this));
}

G09InequalityConstraint::G09InequalityConstraint() : TestVectorFunction(7, 4) {}

G09InequalityConstraint::~G09InequalityConstraint() {}

void G09InequalityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {
  const double x1 = 20.0 * x[0] - 10.0;
  const double x2 = 20.0 * x[1] - 10.0;
  const double x3 = 20.0 * x[2] - 10.0;
  const double x4 = 20.0 * x[3] - 10.0;
  const double x5 = 20.0 * x[4] - 10.0;
  const double x6 = 20.0 * x[5] - 10.0;
  const double x7 = 20.0 * x[6] - 10.0;

  value[0] = -127.0 + 2.0 * x1 * x1 + 3.0 * std::pow(x2, 4.0) + x3 + 4.0 * x4 * x4 + 5.0 * x5;
  value[1] = -282.0 + 7.0 * x1 + 3.0 * x2 + 10.0 * x3 * x3 + x4 - x5;
  value[2] = -196.0 + 23.0 * x1 + x2 * x2 + 6.0 * x6 * x6 - 8.0 * x7;
  value[3] = 4.0 * x1 * x1 + x2 * x2 - 3.0 * x1 * x2 + 2.0 * x3 * x3 + 5.0 * x6 - 11.0 * x7;
}

void G09InequalityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G09InequalityConstraint(*this));
}

G09EqualityConstraint::G09EqualityConstraint() : TestVectorFunction(7, 0) {}

G09EqualityConstraint::~G09EqualityConstraint() {}

void G09EqualityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {}

void G09EqualityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G09EqualityConstraint(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

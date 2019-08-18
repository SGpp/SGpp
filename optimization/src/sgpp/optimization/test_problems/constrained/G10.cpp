// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/G10.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

G10::G10() : ConstrainedTestProblem(8), f(), g(), h() {}

G10::~G10() {}

TestScalarFunction& G10::getObjectiveFunction() { return f; }

TestVectorFunction& G10::getInequalityConstraintFunction() { return g; }

TestVectorFunction& G10::getEqualityConstraintFunction() { return h; }

double G10::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(8);
  x[0] = 0.0484158282828283;
  x[1] = 0.0399936666666667;
  x[2] = 0.456674555555556;
  x[3] = 0.173754949494950;
  x[4] = 0.288483333333333;
  x[5] = 0.210080707070707;
  x[6] = 0.279208282828283;
  x[7] = 0.389492828282828;
  return 7049.3307;
}

G10Objective::G10Objective() : TestScalarFunction(8) {}

G10Objective::~G10Objective() {}

double G10Objective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 9900.0 * x[0] + 100.0;
  const double x2 = 9000.0 * x[1] + 1000.0;
  const double x3 = 9000.0 * x[2] + 1000.0;
  // const double x4 = 990.0 * x[3] + 10.0;
  // const double x5 = 990.0 * x[4] + 10.0;
  // const double x6 = 990.0 * x[5] + 10.0;
  // const double x7 = 990.0 * x[6] + 10.0;
  // const double x8 = 990.0 * x[7] + 10.0;

  return x1 + x2 + x3;
}

void G10Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new G10Objective(*this));
}

G10InequalityConstraint::G10InequalityConstraint() : TestVectorFunction(8, 6) {}

G10InequalityConstraint::~G10InequalityConstraint() {}

void G10InequalityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {
  const double x1 = 9900.0 * x[0] + 100.0;
  const double x2 = 9000.0 * x[1] + 1000.0;
  const double x3 = 9000.0 * x[2] + 1000.0;
  const double x4 = 990.0 * x[3] + 10.0;
  const double x5 = 990.0 * x[4] + 10.0;
  const double x6 = 990.0 * x[5] + 10.0;
  const double x7 = 990.0 * x[6] + 10.0;
  const double x8 = 990.0 * x[7] + 10.0;

  value[0] = -1.0 + (x4 + x6) / 400.0;
  value[1] = -1.0 + (x5 + x7 - x4) / 400.0;
  value[2] = -1.0 + (x8 - x5) / 100.0;
  value[3] = -x1 * x6 + 833.33252 * x4 + 100.0 * x1 - 83333.333;
  value[4] = -x2 * x7 + 1250.0 * x5 + x2 * x4 - 1250.0 * x4;
  value[5] = -x3 * x8 + 1250000.0 + x3 * x5 - 2500.0 * x5;
}

void G10InequalityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G10InequalityConstraint(*this));
}

G10EqualityConstraint::G10EqualityConstraint() : TestVectorFunction(8, 0) {}

G10EqualityConstraint::~G10EqualityConstraint() {}

void G10EqualityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {}

void G10EqualityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G10EqualityConstraint(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

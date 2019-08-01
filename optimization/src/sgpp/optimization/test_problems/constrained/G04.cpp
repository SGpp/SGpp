// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/G04.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

G04::G04() : ConstrainedTestProblem(5), f(), g(), h() {}

G04::~G04() {}

TestScalarFunction& G04::getObjectiveFunction() { return f; }

TestVectorFunction& G04::getInequalityConstraintFunction() { return g; }

TestVectorFunction& G04::getEqualityConstraintFunction() { return h; }

double G04::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(5);
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.166403112537889;
  x[3] = 1.0;
  x[4] = 0.543100716988222;
  return -30665.538671783201;
}

bool G04::isDisplacementFeasible() {
  displacement[0] = 0.0;
  displacement[1] = 0.0;
  displacement[3] = 0.0;
  return ConstrainedTestProblem::isDisplacementFeasible();
}

G04Objective::G04Objective() : TestScalarFunction(5) {}

G04Objective::~G04Objective() {}

double G04Objective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 24.0 * x[0] + 78.0;
  // const double x2 = 12.0 * x[1] + 33.0;
  const double x3 = 18.0 * x[2] + 27.0;
  // const double x4 = 18.0 * x[3] + 27.0;
  const double x5 = 18.0 * x[4] + 27.0;

  return 5.3578547 * x3 * x3 + 0.8356891 * x1 * x5 + 37.293239 * x1 - 40792.141;
}

void G04Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new G04Objective(*this));
}

G04InequalityConstraint::G04InequalityConstraint() : TestVectorFunction(5, 6) {}

G04InequalityConstraint::~G04InequalityConstraint() {}

void G04InequalityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {
  const double x1 = 24.0 * x[0] + 78.0;
  const double x2 = 12.0 * x[1] + 33.0;
  const double x3 = 18.0 * x[2] + 27.0;
  const double x4 = 18.0 * x[3] + 27.0;
  const double x5 = 18.0 * x[4] + 27.0;

  value[0] = 85.334407 + 0.0056858 * x2 * x5 + 0.0006262 * x1 * x4 - 0.0022053 * x3 * x5 - 92.0;
  value[1] = -85.334407 - 0.0056858 * x2 * x5 - 0.0006262 * x1 * x4 + 0.0022053 * x3 * x5;
  value[2] = 80.51249 + 0.0071317 * x2 * x5 + 0.0029955 * x1 * x2 + 0.0021813 * x3 * x3 - 110.0;
  value[3] = -80.51249 - 0.0071317 * x2 * x5 - 0.0029955 * x1 * x2 - 0.0021813 * x3 * x3 + 90.0;
  value[4] = 9.300961 + 0.0047026 * x3 * x5 + 0.0012547 * x1 * x3 + 0.0019085 * x3 * x4 - 25.0;
  value[5] = -9.300961 - 0.0047026 * x3 * x5 - 0.0012547 * x1 * x3 - 0.0019085 * x3 * x4 + 20.0;
}

void G04InequalityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G04InequalityConstraint(*this));
}

G04EqualityConstraint::G04EqualityConstraint() : TestVectorFunction(5, 0) {}

G04EqualityConstraint::~G04EqualityConstraint() {}

void G04EqualityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {}

void G04EqualityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new G04EqualityConstraint(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

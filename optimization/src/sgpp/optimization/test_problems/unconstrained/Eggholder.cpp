// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Eggholder.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Eggholder::Eggholder() : UnconstrainedTestProblem(2), f() {}

Eggholder::~Eggholder() {}

TestScalarFunction& Eggholder::getObjectiveFunction() { return f; }

double Eggholder::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 1.0;
  x[1] = 0.8947577;
  return -959.6406627136321;
}

bool Eggholder::isDisplacementFeasible() {
  displacement[0] = 0.0;
  return UnconstrainedTestProblem::isDisplacementFeasible();
}

EggholderObjective::EggholderObjective() : TestScalarFunction(2) {}

EggholderObjective::~EggholderObjective() {}

double EggholderObjective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 1024.0 * x[0] - 512.0;
  const double x2 = 1024.0 * x[1] - 512.0;

  return -(x2 + 47.0) * std::sin(std::sqrt(std::abs(x1 / 2.0 + x2 + 47.0))) -
         x1 * std::sin(std::sqrt(std::abs(x1 - (x2 + 47.0))));
}

void EggholderObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new EggholderObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

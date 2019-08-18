// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Schwefel06.hpp>
#include <algorithm>
#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Schwefel06::Schwefel06() : UnconstrainedTestProblem(2), f() {}

Schwefel06::~Schwefel06() {}

TestScalarFunction& Schwefel06::getObjectiveFunction() { return f; }

double Schwefel06::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.7;
  x[1] = 0.9;
  return 0.0;
}

Schwefel06Objective::Schwefel06Objective() : TestScalarFunction(2) {}

Schwefel06Objective::~Schwefel06Objective() {}

double Schwefel06Objective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 10.0 * x[0] - 6.0;
  const double x2 = 10.0 * x[1] - 6.0;

  return std::max(std::abs(x1 + 2.0 * x2 - 7), std::abs(2.0 * x1 + x2 - 5));
}

void Schwefel06Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new Schwefel06Objective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

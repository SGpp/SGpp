// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/IncreasingPower.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

IncreasingPower::IncreasingPower(size_t d) : UnconstrainedTestProblem(d), f(d) {}

IncreasingPower::~IncreasingPower() {}

TestScalarFunction& IncreasingPower::getObjectiveFunction() { return f; }

double IncreasingPower::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);

  for (size_t t = 0; t < d; t++) {
    x[t] = 1.0 / std::pow(2.0, t + 1);
  }

  return 0.0;
}

IncreasingPowerObjective::IncreasingPowerObjective(size_t d) : TestScalarFunction(d) {}

IncreasingPowerObjective::~IncreasingPowerObjective() {}

double IncreasingPowerObjective::evalUndisplaced(const base::DataVector& x) {
  double result = 0.0;

  for (size_t t = 0; t < d; t++) {
    result += std::pow(std::abs(x[t] - 1.0 / std::pow(2.0, t + 1)), t + 2);
  }

  return result;
}

void IncreasingPowerObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new IncreasingPowerObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

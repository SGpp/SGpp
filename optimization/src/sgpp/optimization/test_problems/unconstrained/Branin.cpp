// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Branin.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Branin::Branin() : UnconstrainedTestProblem(2), f() {}

Branin::~Branin() {}

TestScalarFunction& Branin::getObjectiveFunction() { return f; }

double Branin::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.5427728435726528825641;
  x[1] = 0.151666666666666666666666667;
  return 5.0 / (4.0 * M_PI);
}

BraninObjective::BraninObjective() : TestScalarFunction(2) {}

BraninObjective::~BraninObjective() {}

double BraninObjective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 15.0 * x[0] - 5.0;
  const double x2 = 15.0 * x[1];
  const double tmp = x2 - 5.1 * x1 * x1 / (4.0 * M_PI * M_PI) + 5.0 * x1 / M_PI - 6.0;

  return tmp * tmp + 10.0 * (1.0 - 1.0 / (8.0 * M_PI)) * std::cos(x1) + 10.0;
}

void BraninObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new BraninObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

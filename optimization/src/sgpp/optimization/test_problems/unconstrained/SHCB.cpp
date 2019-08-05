// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/SHCB.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

SHCB::SHCB() : UnconstrainedTestProblem(2), f() {}

SHCB::~SHCB() {}

TestScalarFunction& SHCB::getObjectiveFunction() { return f; }

double SHCB::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.50898420131003180624224905;
  x[1] = 0.42873435969792603666027341858;
  return -1.031628453489877;
}

SHCBObjective::SHCBObjective() : TestScalarFunction(2) {}

SHCBObjective::~SHCBObjective() {}

double SHCBObjective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 10.0 * x[0] - 5.0;
  const double x2 = 10.0 * x[1] - 5.0;

  return x1 * x1 * (4.0 - 2.1 * x1 * x1 + x1 * x1 * x1 * x1 / 3.0) + x1 * x2 +
         4.0 * x2 * x2 * (x2 * x2 - 1.0);
}

void SHCBObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new SHCBObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

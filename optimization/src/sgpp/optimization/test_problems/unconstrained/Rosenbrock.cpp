// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Rosenbrock.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Rosenbrock::Rosenbrock(size_t d) : UnconstrainedTestProblem(d), f(d) {}

Rosenbrock::~Rosenbrock() {}

TestScalarFunction& Rosenbrock::getObjectiveFunction() { return f; }

double Rosenbrock::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  x.setAll(0.4);
  return 0.0;
}

RosenbrockObjective::RosenbrockObjective(size_t d) : TestScalarFunction(d) {}

RosenbrockObjective::~RosenbrockObjective() {}

double RosenbrockObjective::evalUndisplaced(const base::DataVector& x) {
  double result = 0.0;

  double xt = 15.0 * x[0] - 5.0;

  for (size_t t = 1; t < d; t++) {
    const double xtm1 = xt;
    xt = 15.0 * x[t] - 5.0;

    const double tmp1 = xt - xtm1 * xtm1;
    const double tmp2 = 1.0 - xtm1;
    result += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
  }

  return result;
}

void RosenbrockObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new RosenbrockObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

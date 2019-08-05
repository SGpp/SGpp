// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Rastrigin.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Rastrigin::Rastrigin(size_t d) : UnconstrainedTestProblem(d), f(d) {}

Rastrigin::~Rastrigin() {}

TestScalarFunction& Rastrigin::getObjectiveFunction() { return f; }

double Rastrigin::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  x.setAll(0.2);
  return 0.0;
}

RastriginObjective::RastriginObjective(size_t d) : TestScalarFunction(d) {}

RastriginObjective::~RastriginObjective() {}

double RastriginObjective::evalUndisplaced(const base::DataVector& x) {
  double result = 10.0 * static_cast<double>(d);

  for (size_t t = 0; t < d; t++) {
    const double xt = 10.0 * x[t] - 2.0;
    result += xt * xt - 10.0 * std::cos(2 * M_PI * xt);
  }

  return result;
}

void RastriginObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new RastriginObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

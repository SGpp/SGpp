// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Griewank.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Griewank::Griewank(size_t d) : UnconstrainedTestProblem(d), f(d) {}

Griewank::~Griewank() {}

TestScalarFunction& Griewank::getObjectiveFunction() { return f; }

double Griewank::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  x.setAll(0.5);
  return 0.0;
}

GriewankObjective::GriewankObjective(size_t d) : TestScalarFunction(d) {}

GriewankObjective::~GriewankObjective() {}

double GriewankObjective::evalUndisplaced(const base::DataVector& x) {
  double result = 1.0;
  double tmp = 1.0;

  for (size_t t = 0; t < d; t++) {
    const double xt = 1200.0 * x[t] - 600.0;
    result += xt * xt / 4000.0;
    tmp *= std::cos(xt / std::sqrt(static_cast<double>(t + 1)));
  }

  result -= tmp;
  return result;
}

void GriewankObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new GriewankObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

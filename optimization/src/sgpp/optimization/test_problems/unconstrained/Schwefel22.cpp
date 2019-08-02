// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Schwefel22.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Schwefel22::Schwefel22(size_t d) : UnconstrainedTestProblem(d), f(d) {}

Schwefel22::~Schwefel22() {}

TestScalarFunction& Schwefel22::getObjectiveFunction() { return f; }

double Schwefel22::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  x.setAll(0.3);
  return 0.0;
}

Schwefel22Objective::Schwefel22Objective(size_t d) : TestScalarFunction(d) {}

Schwefel22Objective::~Schwefel22Objective() {}

double Schwefel22Objective::evalUndisplaced(const base::DataVector& x) {
  double result = 0.0;
  double product = 1.0;

  for (size_t t = 0; t < d; t++) {
    const double xt = std::abs(10.0 * x[t] - 3.0);
    result += xt;
    product *= xt;
  }

  result += product;

  return result;
}

void Schwefel22Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new Schwefel22Objective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

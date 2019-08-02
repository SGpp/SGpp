// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Alpine02.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Alpine02::Alpine02(size_t d) : UnconstrainedTestProblem(d), f(d) {}

Alpine02::~Alpine02() {}

TestScalarFunction& Alpine02::getObjectiveFunction() { return f; }

double Alpine02::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  x.setAll(0.739631585583250);
  return -std::pow(std::sqrt(8.0 * x[0] + 2.0) * std::sin(8.0 * x[0] + 2.0),
                   static_cast<double>(d));
}

Alpine02Objective::Alpine02Objective(size_t d) : TestScalarFunction(d) {}

Alpine02Objective::~Alpine02Objective() {}

double Alpine02Objective::evalUndisplaced(const base::DataVector& x) {
  double result = -1.0;

  for (size_t t = 0; t < d; t++) {
    const double xt = 8.0 * x[t] + 2.0;
    result *= std::sqrt(xt) * std::sin(xt);
  }

  return result;
}

void Alpine02Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new Alpine02Objective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

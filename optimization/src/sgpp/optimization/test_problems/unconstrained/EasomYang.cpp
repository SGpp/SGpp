// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/EasomYang.hpp>
#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

EasomYang::EasomYang(size_t d) : UnconstrainedTestProblem(d), f(d) {}

EasomYang::~EasomYang() {}

TestScalarFunction& EasomYang::getObjectiveFunction() { return f; }

double EasomYang::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  x.setAll(0.75);
  return -1.0;
}

EasomYangObjective::EasomYangObjective(size_t d) : TestScalarFunction(d) {}

EasomYangObjective::~EasomYangObjective() {}

double EasomYangObjective::evalUndisplaced(const base::DataVector& x) {
  double sum = 0.0;
  double product = 1.0;

  for (size_t t = 0; t < d; t++) {
    const double xt = 2.0 * M_PI * (2.0 * x[t] - 1.0);
    const double diff = xt - M_PI;
    sum += diff * diff;
    product *= -std::cos(xt);
  }

  return -std::exp(-sum) * product;
}

void EasomYangObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new EasomYangObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

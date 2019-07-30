// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/BubbleWrap.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

BubbleWrap::BubbleWrap(size_t d) : UnconstrainedTestProblem(d), f(d) {}

BubbleWrap::~BubbleWrap() {}

TestScalarFunction& BubbleWrap::getObjectiveFunction() { return f; }

double BubbleWrap::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  x.setAll(1.0 / 3.0);
  return 0.0;
}

BubbleWrapObjective::BubbleWrapObjective(size_t d) : TestScalarFunction(d) {}

BubbleWrapObjective::~BubbleWrapObjective() {}

double BubbleWrapObjective::evalUndisplaced(const base::DataVector& x) {
  double product = 1.0;

  for (size_t t = 0; t < d; t++) {
    const double xt = 1.5 * x[t] - 0.5;
    product *= 0.9 - std::abs(xt) + 0.1 * std::cos(10.0 * M_PI * xt);
  }

  return 1.0 - product;
}

void BubbleWrapObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new BubbleWrapObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

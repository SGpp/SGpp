// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/GoldsteinPrice.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

GoldsteinPrice::GoldsteinPrice() : UnconstrainedTestProblem(2), f() {}

GoldsteinPrice::~GoldsteinPrice() {}

TestScalarFunction& GoldsteinPrice::getObjectiveFunction() { return f; }

double GoldsteinPrice::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.5;
  x[1] = 0.25;
  return 0.0003;
}

GoldsteinPriceObjective::GoldsteinPriceObjective() : TestScalarFunction(2) {}

GoldsteinPriceObjective::~GoldsteinPriceObjective() {}

double GoldsteinPriceObjective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 4.0 * x[0] - 2.0;
  const double x2 = 4.0 * x[1] - 2.0;

  return (1.0 +
          (x1 + x2 + 1.0) * (x1 + x2 + 1.0) *
              (19.0 - 14.0 * x1 + 3.0 * x1 * x1 - 14.0 * x2 + 6.0 * x1 * x2 +
                  3.0 * x2 * x2)) *
         (30.0 +
          (2.0 * x1 - 3.0 * x2) * (2.0 * x1 - 3.0 * x2) *
              (18.0 - 32.0 * x1 + 12.0 * x1 * x1 + 48.0 * x2 - 36.0 * x1 * x2 +
                  27.0 * x2 * x2)) * 1e-4;
}

void GoldsteinPriceObjective::clone(std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::unique_ptr<ScalarFunction>(

      new GoldsteinPriceObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Himmelblau.hpp>

#include <cmath>

namespace SGPP {
namespace optimization {
namespace test_problems {

Himmelblau::Himmelblau() : UnconstrainedTestProblem(2), f() {}

Himmelblau::~Himmelblau() {}

TestScalarFunction& Himmelblau::getObjectiveFunction() { return f; }

float_t Himmelblau::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.8;
  x[1] = 0.7;
  return 0.0;
}

HimmelblauObjective::HimmelblauObjective() : TestScalarFunction(2) {}

HimmelblauObjective::~HimmelblauObjective() {}

float_t HimmelblauObjective::evalUndisplaced(const base::DataVector& x) {
  const float_t x1 = 10.0 * x[0] - 5.0;
  const float_t x2 = 10.0 * x[1] - 5.0;

  return (x1 * x1 + x2 - 11.0) * (x1 * x1 + x2 - 11.0) +
         (x1 + x2 * x2 - 7.0) * (x1 + x2 * x2 - 7.0);
}

void HimmelblauObjective::clone(std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::unique_ptr<ScalarFunction>(new HimmelblauObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace SGPP

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Hartman3.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Hartman3::Hartman3() : UnconstrainedTestProblem(3), f() {}

Hartman3::~Hartman3() {}

TestScalarFunction& Hartman3::getObjectiveFunction() { return f; }

double Hartman3::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(3);
  x[0] = 0.1146398;
  x[1] = 0.5556488;
  x[2] = 0.8525470;
  return -3.862784507551574;
}

Hartman3Objective::Hartman3Objective() : TestScalarFunction(3) {}

Hartman3Objective::~Hartman3Objective() {}

double Hartman3Objective::evalUndisplaced(const base::DataVector& x) {
  return -1.0 * std::exp(-3.0 * (x[0] - 0.3689) * (x[0] - 0.3689) -
                         10.0 * (x[1] - 0.1170) * (x[1] - 0.1170) -
                         30.0 * (x[2] - 0.2673) * (x[2] - 0.2673)) -
         1.2 * std::exp(-0.1 * (x[0] - 0.4699) * (x[0] - 0.4699) -
                        10.0 * (x[1] - 0.4387) * (x[1] - 0.4387) -
                        35.0 * (x[2] - 0.7470) * (x[2] - 0.7470)) -
         3.0 * std::exp(-3.0 * (x[0] - 0.1091) * (x[0] - 0.1091) -
                        10.0 * (x[1] - 0.8732) * (x[1] - 0.8732) -
                        30.0 * (x[2] - 0.5547) * (x[2] - 0.5547)) -
         3.2 * std::exp(-0.1 * (x[0] - 0.0382) * (x[0] - 0.0382) -
                        10.0 * (x[1] - 0.5743) * (x[1] - 0.5743) -
                        35.0 * (x[2] - 0.8828) * (x[2] - 0.8828));
}

void Hartman3Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new Hartman3Objective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

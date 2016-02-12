// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/HoelderTable.hpp>

#include <cmath>

namespace SGPP {
namespace optimization {
namespace test_problems {

HoelderTable::HoelderTable() :
  UnconstrainedTestProblem(2),
  f() {
}

HoelderTable::~HoelderTable() {
}

TestScalarFunction& HoelderTable::getObjectiveFunction() {
  return f;
}

float_t HoelderTable::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.90275115;
  x[1] = 0.9832295;
  return -19.208502567884491;
}

bool HoelderTable::isDisplacementFeasible() {
  if ((displacement[0] > 0.005) || (displacement[0] < -0.005) ||
      (displacement[1] > 0.01) || (displacement[1] < -0.01)) {
    return false;
  }

  return UnconstrainedTestProblem::isDisplacementFeasible();
}

HoelderTableObjective::HoelderTableObjective() :
  TestScalarFunction(2) {
}

HoelderTableObjective::~HoelderTableObjective() {
}

float_t HoelderTableObjective::evalUndisplaced(
  const base::DataVector& x) {
  const float_t x1 = 20.0 * x[0] - 10.0;
  const float_t x2 = 20.0 * x[1] - 10.0;

  return -std::abs(
           std::sin(x1) * std::cos(x2) * std::exp(
             std::abs(1.0 - std::sqrt(x1 * x1 + x2 * x2) / M_PI)));
}

void HoelderTableObjective::clone(
  std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::unique_ptr<ScalarFunction>(

            new HoelderTableObjective(*this));
}

}
}
}

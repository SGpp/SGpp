// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Michalewicz.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Michalewicz::Michalewicz() : UnconstrainedTestProblem(2), f() {}

Michalewicz::~Michalewicz() {}

TestScalarFunction& Michalewicz::getObjectiveFunction() { return f; }

double Michalewicz::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.440581104135123915;
  x[1] = M_PI / 10.0;
  return -1.801303410098554;
}

MichalewiczObjective::MichalewiczObjective() : TestScalarFunction(2) {}

MichalewiczObjective::~MichalewiczObjective() {}

double MichalewiczObjective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 5.0 * x[0];
  const double x2 = 5.0 * x[1];

  return -std::sin(x1) * std::pow(std::sin(x1 * x1 / M_PI), 20.0) -
         std::sin(x2) * std::pow(std::sin(2.0 * x2 * x2 / M_PI), 20.0);
}

void MichalewiczObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new MichalewiczObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

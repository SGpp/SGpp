// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Beale.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Beale::Beale() : UnconstrainedTestProblem(2), f() {}

Beale::~Beale() {}

TestScalarFunction& Beale::getObjectiveFunction() { return f; }

double Beale::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.8;
  x[1] = 0.55;
  return 0.0;
}

BealeObjective::BealeObjective() : TestScalarFunction(2) {}

BealeObjective::~BealeObjective() {}

double BealeObjective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 10.0 * x[0] - 5.0;
  const double x2 = 10.0 * x[1] - 5.0;
  const double tmp1 = 1.5 - x1 * (1.0 - x2);
  const double tmp2 = 2.25 - x1 * (1.0 - x2 * x2);
  const double tmp3 = 2.625 - x1 * (1.0 - x2 * x2 * x2);

  return tmp1 * tmp1 + tmp2 * tmp2 + tmp3 * tmp3;
}

void BealeObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new BealeObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

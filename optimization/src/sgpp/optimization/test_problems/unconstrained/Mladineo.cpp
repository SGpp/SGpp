// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Mladineo.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Mladineo::Mladineo() : UnconstrainedTestProblem(2), f() {}

Mladineo::~Mladineo() {}

TestScalarFunction& Mladineo::getObjectiveFunction() { return f; }

double Mladineo::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.001542464646465;
  x[1] = 0.004449101010101;
  return 1.701830178238994e-04;
}

bool Mladineo::isDisplacementFeasible() {
  if ((displacement[0] > 0) || (displacement[0] < -0.01) || (displacement[1] > 0) ||
      (displacement[1] < -0.01)) {
    return false;
  }

  return UnconstrainedTestProblem::isDisplacementFeasible();
}

MladineoObjective::MladineoObjective() : TestScalarFunction(2) {}

MladineoObjective::~MladineoObjective() {}

double MladineoObjective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 0.99 * x[0] + 0.01;
  const double x2 = 0.99 * x[1] + 0.01;

  return 1.0 + (x1 * x1 + x2 * x2) / 2.0 -
         std::cos(10.0 * std::log(2.0 * x1)) * std::cos(10.0 * std::log(3.0 * x2));
}

void MladineoObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new MladineoObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

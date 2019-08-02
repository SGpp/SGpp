// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Branin02.hpp>
#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Branin02::Branin02() : UnconstrainedTestProblem(2), f() {}

Branin02::~Branin02() {}

TestScalarFunction& Branin02::getObjectiveFunction() { return f; }

double Branin02::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.0901505787598000;
  x[1] = 0.876312894266000;
  return 5.55891440389382;
}

Branin02Objective::Branin02Objective() : TestScalarFunction(2) {}

Branin02Objective::~Branin02Objective() {}

double Branin02Objective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 20.0 * x[0] - 5.0;
  const double x2 = 20.0 * x[1] - 5.0;
  const double tmp = x2 - 5.1 * x1 * x1 / (4.0 * M_PI * M_PI) + 5.0 * x1 / M_PI - 6.0;

  return tmp * tmp + 10.0 * (1.0 - 1.0 / (8.0 * M_PI)) * std::cos(x1) * std::cos(x2) +
      std::log(x1 * x1 + x2 * x2 + 1.0) + 10.0;
}

void Branin02Objective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new Branin02Objective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

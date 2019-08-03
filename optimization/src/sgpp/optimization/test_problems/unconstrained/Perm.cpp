// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Perm.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Perm::Perm(size_t d) : UnconstrainedTestProblem(d), f(d) {}

Perm::~Perm() {}

TestScalarFunction& Perm::getObjectiveFunction() { return f; }

double Perm::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  const double dDbl = static_cast<double>(d);

  for (size_t t = 0; t < d; t++) {
    x[t] = 0.5 * (static_cast<double>(t + 1) / dDbl + 1.0);
  }

  return 0.0;
}

PermObjective::PermObjective(size_t d) : TestScalarFunction(d) {}

PermObjective::~PermObjective() {}

double PermObjective::evalUndisplaced(const base::DataVector& x) {
  double result = 0.0;
  const double dDbl = static_cast<double>(d);

  for (size_t i = 0; i < d; i++) {
    const double iDbl = static_cast<double>(i + 1);
    double innerSum = 0.0;

    for (size_t t = 0; t < d; t++) {
      const double xt = dDbl * (2.0 * x[t] - 1.0);
      const double tDbl = static_cast<double>(t + 1);

      innerSum += (std::pow(tDbl, iDbl) + 1.0) * (std::pow(xt / tDbl, iDbl) - 1.0);
    }

    result += std::pow(innerSum, 2.0);
  }

  return result;
}

void PermObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new PermObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

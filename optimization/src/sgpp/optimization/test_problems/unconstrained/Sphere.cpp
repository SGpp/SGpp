// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Sphere.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Sphere::Sphere(size_t d) : UnconstrainedTestProblem(d), f(d) {}

Sphere::~Sphere() {}

TestScalarFunction& Sphere::getObjectiveFunction() { return f; }

double Sphere::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  x.setAll(0.1);
  return 0.0;
}

SphereObjective::SphereObjective(size_t d) : TestScalarFunction(d) {}

SphereObjective::~SphereObjective() {}

double SphereObjective::evalUndisplaced(const base::DataVector& x) {
  double result = 0.0;

  for (size_t t = 0; t < d; t++) {
    const double xt = 10.0 * x[t] - 1.0;
    result += xt * xt;
  }

  return result;
}

void SphereObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new SphereObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

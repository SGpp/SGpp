// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Sphere.hpp>

#include <cmath>

namespace SGPP {
namespace optimization {
namespace test_problems {

Sphere::Sphere(size_t d) : UnconstrainedTestProblem(d), f(d) {}

Sphere::~Sphere() {}

TestScalarFunction& Sphere::getObjectiveFunction() { return f; }

float_t Sphere::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  x.setAll(0.1);
  return 0.0;
}

SphereObjective::SphereObjective(size_t d) : TestScalarFunction(d) {}

SphereObjective::~SphereObjective() {}

float_t SphereObjective::evalUndisplaced(const base::DataVector& x) {
  float_t result = 0.0;

  for (size_t t = 0; t < d; t++) {
    const float_t xt = 10.0 * x[t] - 1.0;
    result += xt * xt;
  }

  return result;
}

void SphereObjective::clone(std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::unique_ptr<ScalarFunction>(new SphereObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace SGPP

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/Soland.hpp>
#include <cmath>
#include "../../../../../../base/src/sgpp/base/function/vector/EmptyVectorFunction.hpp"

namespace sgpp {
namespace optimization {
namespace test_problems {

Soland::Soland() : ConstrainedTestProblem(2), f(), g(), h() {}

Soland::~Soland() {}

TestScalarFunction& Soland::getObjectiveFunction() { return f; }

TestVectorFunction& Soland::getInequalityConstraintFunction() { return g; }

TestVectorFunction& Soland::getEqualityConstraintFunction() { return h; }

double Soland::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.358768096599359;
  x[1] = 0.489947363788600;
  return -16.738893184394644;
}

SolandObjective::SolandObjective() : TestScalarFunction(2) {}

SolandObjective::~SolandObjective() {}

double SolandObjective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 2.0 * x[0];
  const double x2 = 3.0 * x[1];

  return -12.0 * x1 - 7.0 * x2 + x2 * x2;
}

void SolandObjective::clone(std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::unique_ptr<ScalarFunction>(new SolandObjective(*this));
}

SolandInequalityConstraint::SolandInequalityConstraint() : TestVectorFunction(2, 0) {}

SolandInequalityConstraint::~SolandInequalityConstraint() {}

void SolandInequalityConstraint::evalUndisplaced(const base::DataVector& x,
                                                 base::DataVector& value) {}

void SolandInequalityConstraint::clone(std::unique_ptr<VectorFunction>& clone) const {
  clone = std::unique_ptr<VectorFunction>(new SolandInequalityConstraint(*this));
}

SolandEqualityConstraint::SolandEqualityConstraint() : TestVectorFunction(2, 1) {}

SolandEqualityConstraint::~SolandEqualityConstraint() {}

void SolandEqualityConstraint::evalUndisplaced(const base::DataVector& x, base::DataVector& value) {
  const double x1 = 2.0 * x[0];
  const double x2 = 3.0 * x[1];

  value[0] = -2.0 * std::pow(x1, 4.0) + 2.0 - x2;
}

void SolandEqualityConstraint::clone(std::unique_ptr<VectorFunction>& clone) const {
  clone = std::unique_ptr<VectorFunction>(new SolandEqualityConstraint(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

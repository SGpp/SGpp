// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/function/vector/EmptyVectorFunction.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/constrained/Floudas.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

Floudas::Floudas() : ConstrainedTestProblem(2), f(), g(), h() {}

Floudas::~Floudas() {}

TestScalarFunction& Floudas::getObjectiveFunction() { return f; }

TestVectorFunction& Floudas::getInequalityConstraintFunction() { return g; }

TestVectorFunction& Floudas::getEqualityConstraintFunction() { return h; }

double Floudas::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(2);
  x[0] = 0.776506732492537;
  x[1] = 0.794623268529425;
  return -5.50801327159531;
}

FloudasObjective::FloudasObjective() : TestScalarFunction(2) {}

FloudasObjective::~FloudasObjective() {}

double FloudasObjective::evalUndisplaced(const base::DataVector& x) {
  const double x1 = 3.0 * x[0];
  const double x2 = 4.0 * x[1];

  return -x1 - x2;
}

void FloudasObjective::clone(std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::unique_ptr<ScalarFunction>(new FloudasObjective(*this));
}

FloudasInequalityConstraint::FloudasInequalityConstraint() : TestVectorFunction(2, 2) {}

FloudasInequalityConstraint::~FloudasInequalityConstraint() {}

void FloudasInequalityConstraint::evalUndisplaced(const base::DataVector& x,
                                                  base::DataVector& value) {
  const double x1 = 3.0 * x[0];
  const double x2 = 4.0 * x[1];

  const double x1p2 = x1 * x1;
  const double x1p3 = x1p2 * x1;
  const double x1p4 = x1p3 * x1;

  value[0] = x2 - (2.0 * x1p4 - 8.0 * x1p3 + 8.0 * x1p2 + 2.0);
  value[1] = x2 - (4.0 * x1p4 - 32.0 * x1p3 + 88.0 * x1p2 - 96.0 * x1 + 36.0);
}

void FloudasInequalityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new FloudasInequalityConstraint(*this));
}

FloudasEqualityConstraint::FloudasEqualityConstraint() : TestVectorFunction(2, 0) {}

FloudasEqualityConstraint::~FloudasEqualityConstraint() {}

void FloudasEqualityConstraint::evalUndisplaced(const base::DataVector& x,
                                                base::DataVector& value) {}

void FloudasEqualityConstraint::clone(std::unique_ptr<base::VectorFunction>& clone) const {
  clone = std::unique_ptr<base::VectorFunction>(new FloudasEqualityConstraint(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Ackley.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

Ackley::Ackley(size_t d) : UnconstrainedTestProblem(d), f(d) {}

Ackley::~Ackley() {}

TestScalarFunction& Ackley::getObjectiveFunction() { return f; }

double Ackley::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  x.setAll(0.0948903972968234);
  return 6.559645375627878;
}

AckleyObjective::AckleyObjective(size_t d) : TestScalarFunction(d) {}

AckleyObjective::~AckleyObjective() {}

double AckleyObjective::evalUndisplaced(const base::DataVector& x) {
  double result = 0.0;

  double arg1 = 0.0;
  double arg2 = 0.0;

  for (size_t t = 0; t < d; t++) {
    const double xt = 5.0 * x[t] + 1.5;
    arg1 += xt * xt;
    arg2 += std::cos(2.0 * M_PI * xt);
  }

  result = 20.0 * (1.0 - std::exp(-0.2 * std::sqrt(arg1 / static_cast<double>(d))));
  result += M_E - std::exp(arg2 / static_cast<double>(d));

  return result;
}

void AckleyObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new AckleyObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

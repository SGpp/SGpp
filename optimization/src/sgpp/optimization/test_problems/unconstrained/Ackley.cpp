// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Ackley.hpp>

#include <cmath>

namespace SGPP {
namespace optimization {
namespace test_problems {

Ackley::Ackley(size_t d) :
  UnconstrainedTestProblem(d),
  f(d) {
}

Ackley::~Ackley() {
}

TestScalarFunction& Ackley::getObjectiveFunction() {
  return f;
}

float_t Ackley::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  x.setAll(0.1);
  return 0.0;
}

AckleyObjective::AckleyObjective(size_t d) :
  TestScalarFunction(d) {
}

AckleyObjective::~AckleyObjective() {
}

float_t AckleyObjective::evalUndisplaced(
  const base::DataVector& x) {
  float_t result = 0.0;

  float_t arg1 = 0.0;
  float_t arg2 = 0.0;

  for (size_t t = 0; t < d; t++) {
    const float_t xt = 10.0 * x[t] - 1.0;
    arg1 += xt * xt;
    arg2 += std::cos(2.0 * M_PI * xt);
  }

  result = 20.0 *
           (1.0 -
            std::exp(-0.2 *
                     std::sqrt(arg1 / static_cast<float_t>(d))));
  result += M_E - std::exp(arg2 / static_cast<float_t>(d));

  return result;
}

void AckleyObjective::clone(
  std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::unique_ptr<ScalarFunction>(
            new AckleyObjective(*this));
}

}
}
}

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/TremblingParabola.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

TremblingParabola::TremblingParabola(size_t d, size_t p) : UnconstrainedTestProblem(d), f(d, p) {}

TremblingParabola::~TremblingParabola() {}

TestScalarFunction& TremblingParabola::getObjectiveFunction() { return f; }

double TremblingParabola::getOptimalPointUndisplaced(base::DataVector& x) {
  x.resize(d);
  x.setAll(0.2);
  return 0.0;
}

TremblingParabolaObjective::TremblingParabolaObjective(size_t d, size_t p)
    : TestScalarFunction(d), p(p), bSplineBasis(p), g0(splineTrembling(0.0)) {}

TremblingParabolaObjective::~TremblingParabolaObjective() {}

double TremblingParabolaObjective::evalUndisplaced(const base::DataVector& x) {
  double result = 0.0;

  /*for (size_t t = 0; t < d; t++) {
    const double xt = 50.0 * x[t] - 10.0;
    const double tmp = std::copysign(xt, std::sin(M_PI * xt));
    result += xt * xt / 10.0 + std::abs(xt) * (tmp - std::floor(tmp));
  }*/

  for (size_t t = 0; t < d; t++) {
    const double xt = 20.0 * x[t] - 4.0;
    const double gxt = splineTrembling(xt);
    result += 0.1 * xt * xt + std::abs(xt) * (gxt / g0 + 1.0) * 0.5;
  }

  return result;
}

inline double TremblingParabolaObjective::splineTrembling(double x) const {
  const double pMinus1Halved = (static_cast<double>(p) - 1.0) / 2.0;
  const int kMin = static_cast<int>(std::floor(x - pMinus1Halved));
  const int kMax = static_cast<int>(std::ceil(x + pMinus1Halved));
  double xCur = x - kMin + (pMinus1Halved + 1.0);
  double sign = (kMin % 2 == 0) ? 1.0 : -1.0;
  double result = 0.0;

  for (int k = kMin; k <= kMax; k++) {
    result += sign * bSplineBasis.uniformBSpline(xCur, p);
    xCur -= 1.0;
    sign *= -1.0;
  }

  return result;
}

void TremblingParabolaObjective::clone(std::unique_ptr<base::ScalarFunction>& clone) const {
  clone = std::unique_ptr<base::ScalarFunction>(new TremblingParabolaObjective(*this));
}
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

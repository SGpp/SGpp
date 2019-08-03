// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/LineSearchArmijo.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NLCG.hpp>

#include <limits>
#include <numeric>

namespace sgpp {
namespace optimization {
namespace optimizer {

NLCG::NLCG(const base::ScalarFunction& f, const base::ScalarFunctionGradient& fGradient,
           size_t maxItCount, double beta, double gamma, double tolerance, double epsilon,
           double restartThreshold)
    : UnconstrainedOptimizer(f, &fGradient, nullptr, maxItCount),
      beta(beta),
      gamma(gamma),
      tol(tolerance),
      eps(epsilon),
      alpha(restartThreshold) {
}

NLCG::NLCG(const NLCG& other)
    : UnconstrainedOptimizer(other),
      beta(other.beta),
      gamma(other.gamma),
      tol(other.tol),
      eps(other.eps),
      alpha(other.alpha) {
}

NLCG::~NLCG() {}

void NLCG::optimize() {
  base::Printer::getInstance().printStatusBegin("Optimizing (NLCG)...");

  const size_t d = f->getNumberOfParameters();

  xOpt.resize(0);
  fOpt = std::numeric_limits<double>::quiet_NaN();
  xHist.resize(0, d);
  fHist.resize(0);

  base::DataVector x(x0);
  double fx;
  double fy;

  base::DataVector gradFx(d);
  base::DataVector gradFy(d);
  base::DataVector s(d);
  base::DataVector sNormalized(d);
  base::DataVector y(d);

  fx = fGradient->eval(x0, gradFx);
  double gradFxNorm = gradFx.l2Norm();

  xHist.appendRow(x);
  fHist.append(fx);

  size_t k = 1;

  const double BINDING_TOLERANCE = 1e-6;

  // negated gradient as starting search direction
  for (size_t t = 0; t < d; t++) {
    s[t] = -gradFx[t];
  }

  while (k < N) {
    // exit if norm small enough
    if (gradFxNorm < tol) {
      break;
    }

    for (size_t t = 0; t < d; t++) {
      // is constraint binding?
      // (i.e., we are at the boundary and the search direction points outwards)
      if (((x[t] < BINDING_TOLERANCE) && (s[t] < 0.0)) ||
          ((x[t] > 1.0 - BINDING_TOLERANCE) && (s[t] > 0.0))) {
        // discard variable by setting direction to zero
        s[t] = 0.0;
        gradFx[t] = 0.0;
      }
    }

    // recalculate norm of gradient
    gradFxNorm = gradFx.l2Norm();

    // normalize search direction
    const double sNorm = s.l2Norm();

    for (size_t t = 0; t < d; t++) {
      sNormalized[t] = s[t] / sNorm;
    }

    // line search
    if (!lineSearchArmijo(*f, beta, gamma, tol, eps, x, fx, gradFx, sNormalized, y, k)) {
      // line search failed ==> exit
      // (either a "real" error occurred or the improvement achieved is
      // too small)
      break;
    }

    // calculate gradient and norm
    fy = fGradient->eval(y, gradFy);
    k++;

    for (size_t t = 0; t < d; t++) {
      // is constraint binding?
      // (i.e., we are at the boundary and the negated gradient points outwards)
      if (((y[t] < BINDING_TOLERANCE) && (-gradFy[t] < 0.0)) ||
          ((y[t] > 1.0 - BINDING_TOLERANCE) && (-gradFy[t] > 0.0))) {
        // discard variable by setting gradient to zero
        gradFy[t] = 0.0;
      }
    }

    const double gradFyNorm = gradFy.l2Norm();

    double beta = 0.0;

    // the method is restarted (beta = 0), if the following criterion
    // is *not* met
    if (std::abs(gradFy.dotProduct(gradFx)) / (gradFyNorm * gradFyNorm) < alpha) {
      // Polak-Ribiere coefficient
      for (size_t t = 0; t < d; t++) {
        beta += gradFy[t] * (gradFy[t] - gradFx[t]);
      }

      beta /= gradFxNorm * gradFxNorm;
    }

    // new search direction
    for (size_t t = 0; t < d; t++) {
      s[t] = beta * s[t] - gradFy[t];
    }

    x = y;
    fx = fy;
    gradFx = gradFy;
    gradFxNorm = gradFyNorm;
    xHist.appendRow(x);
    fHist.append(fx);

    // status printing
    base::Printer::getInstance().printStatusUpdate(
        std::to_string(k) + " evaluations, x = " + x.toString() + ", f(x) = " + std::to_string(fx));
  }

  xOpt.resize(d);
  xOpt = x;
  fOpt = fx;
  base::Printer::getInstance().printStatusEnd();
}

double NLCG::getBeta() const { return beta; }

void NLCG::setBeta(double beta) { this->beta = beta; }

double NLCG::getGamma() const { return gamma; }

void NLCG::setGamma(double gamma) { this->gamma = gamma; }

double NLCG::getTolerance() const { return tol; }

void NLCG::setTolerance(double tolerance) { tol = tolerance; }

double NLCG::getEpsilon() const { return eps; }

void NLCG::setEpsilon(double epsilon) { eps = epsilon; }

double NLCG::getRestartThreshold() const { return alpha; }

void NLCG::setRestartThreshold(double restartThreshold) { alpha = restartThreshold; }

void NLCG::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new NLCG(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

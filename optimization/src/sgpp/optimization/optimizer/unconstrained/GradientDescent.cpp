// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/GradientDescent.hpp>
#include <sgpp/optimization/optimizer/unconstrained/LineSearchArmijo.hpp>

namespace sgpp {
namespace optimization {
namespace optimizer {

GradientDescent::GradientDescent(const base::ScalarFunction& f,
                                 const base::ScalarFunctionGradient& fGradient, size_t maxItCount,
                                 double beta, double gamma, double tolerance, double epsilon)
    : UnconstrainedOptimizer(f, maxItCount),
      beta(beta),
      gamma(gamma),
      tol(tolerance),
      eps(epsilon) {
  fGradient.clone(this->fGradient);
}
GradientDescent::GradientDescent(const GradientDescent& other)
    : UnconstrainedOptimizer(other),
      beta(other.beta),
      gamma(other.gamma),
      tol(other.tol),
      eps(other.eps) {
  other.fGradient->clone(fGradient);
}

GradientDescent::~GradientDescent() {}

void GradientDescent::optimize() {
  base::Printer::getInstance().printStatusBegin("Optimizing (gradient descent)...");

  const size_t d = f->getNumberOfParameters();

  xOpt.resize(0);
  fOpt = NAN;
  xHist.resize(0, d);
  fHist.resize(0);

  base::DataVector x(x0);
  double fx = NAN;

  base::DataVector gradFx(d);
  base::DataVector s(d);
  base::DataVector y(d);
  size_t k = 0;

  while (k < N) {
    // calculate gradient and norm
    fx = fGradient->eval(x, gradFx);
    k++;
    double gradFxNorm = gradFx.l2Norm();

    if (k == 1) {
      xHist.appendRow(x);
      fHist.append(fx);
    }

    // exit if norm small enough
    if (gradFxNorm < tol) {
      break;
    }

    // search direction is the normalized negated gradient
    for (size_t t = 0; t < d; t++) {
      s[t] = -gradFx[t] / gradFxNorm;
    }

    // status printing
    base::Printer::getInstance().printStatusUpdate(
        std::to_string(k) + " evaluations, x = " + x.toString() + ", f(x) = " + std::to_string(fx));

    // line search
    if (!lineSearchArmijo(*f, beta, gamma, tol, eps, x, fx, gradFx, s, y, k)) {
      // line search failed ==> exit
      // (either a "real" error occured or the improvement
      // achieved is too small)
      break;
    }

    x = y;
    xHist.appendRow(x);
    fHist.append(fx);
  }

  xOpt.resize(d);
  xOpt = x;
  fOpt = fx;
  base::Printer::getInstance().printStatusEnd();
}

base::ScalarFunctionGradient& GradientDescent::getObjectiveGradient() const { return *fGradient; }

double GradientDescent::getBeta() const { return beta; }

void GradientDescent::setBeta(double beta) { this->beta = beta; }

double GradientDescent::getGamma() const { return gamma; }

void GradientDescent::setGamma(double gamma) { this->gamma = gamma; }

double GradientDescent::getTolerance() const { return tol; }

void GradientDescent::setTolerance(double tolerance) { tol = tolerance; }

double GradientDescent::getEpsilon() const { return eps; }

void GradientDescent::setEpsilon(double epsilon) { eps = epsilon; }

void GradientDescent::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new GradientDescent(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>

namespace sgpp {
namespace optimization {
namespace optimizer {

AdaptiveGradientDescent::AdaptiveGradientDescent(
    const ScalarFunction& f, const ScalarFunctionGradient& fGradient,
    size_t maxItCount, double tolerance,
    double stepSizeIncreaseFactor, double stepSizeDecreaseFactor, double lineSearchAccuracy)
    : UnconstrainedOptimizer(f, maxItCount),
      theta(tolerance),
      rhoAlphaPlus(stepSizeIncreaseFactor),
      rhoAlphaMinus(stepSizeDecreaseFactor),
      rhoLs(lineSearchAccuracy) {
  fGradient.clone(this->fGradient);
}

AdaptiveGradientDescent::AdaptiveGradientDescent(const AdaptiveGradientDescent& other) :
    UnconstrainedOptimizer(other),
    theta(other.theta),
    rhoAlphaPlus(other.rhoAlphaPlus),
    rhoAlphaMinus(other.rhoAlphaMinus),
    rhoLs(other.rhoLs) {
}

AdaptiveGradientDescent::~AdaptiveGradientDescent() {}

void AdaptiveGradientDescent::optimize() {
  Printer::getInstance().printStatusBegin("Optimizing (adaptive gradient descent)...");

  const size_t d = f->getNumberOfParameters();

  xOpt.resize(0);
  fOpt = NAN;
  xHist.resize(0, d);
  fHist.resize(0);

  base::DataVector x(x0);
  double fx = f->eval(x);
  base::DataVector gradFx(d);

  xHist.appendRow(x);
  fHist.append(fx);

  base::DataVector xNew(d);
  double fxNew;

  size_t k = 0;
  double alpha = 1.0;
  base::DataVector dir(d);
  bool inDomain;

  size_t breakIterationCounter = 0;
  const size_t BREAK_ITERATION_COUNTER_MAX = 10;

  while (k < N) {
    // calculate gradient and norm
    fx = fGradient->eval(x, gradFx);
    k++;

    const double gradFxNorm = gradFx.l2Norm();

    if (gradFxNorm == 0.0) {
      break;
    }

    inDomain = true;

    for (size_t t = 0; t < d; t++) {
      // search direction (normalized negated gradient)
      dir[t] = -gradFx[t] / gradFxNorm;
      // new point
      xNew[t] = x[t] + alpha * dir[t];

      if ((xNew[t] < 0.0) || (xNew[t] > 1.0)) {
        inDomain = false;
      }
    }

    // evaluate at new point
    fxNew = (inDomain ? f->eval(xNew) : INFINITY);
    k++;

    // inner product of gradient and search direction
    const double gradFxTimesDir = -gradFxNorm;

    // line search
    while ((fxNew > fx + rhoLs * alpha * gradFxTimesDir) && (alpha > 0.0)) {
      alpha *= rhoAlphaMinus;
      inDomain = true;

      // recalculate new point
      for (size_t t = 0; t < d; t++) {
        xNew[t] = x[t] + alpha * dir[t];

        if ((xNew[t] < 0.0) || (xNew[t] > 1.0)) {
          inDomain = false;
          break;
        }
      }

      // evaluate at new point
      fxNew = (inDomain ? f->eval(xNew) : INFINITY);
      k++;
    }

    // after too many line search steps, alpha will be numerically zero
    if (alpha == 0.0) {
      break;
    }

    // save new point
    x = xNew;
    fx = fxNew;
    xHist.appendRow(x);
    fHist.append(fx);

    // increase step size
    alpha *= rhoAlphaPlus;

    // status printing
    Printer::getInstance().printStatusUpdate(std::to_string(k) + " evaluations, x = " +
                                             x.toString() + ", f(x) = " + std::to_string(fx));

    // stopping criterion:
    // stop if alpha is smaller than tolerance theta
    // in BREAK_ITERATION_COUNTER_MAX consecutive iterations
    if (alpha < theta) {
      breakIterationCounter++;

      if (breakIterationCounter >= BREAK_ITERATION_COUNTER_MAX) {
        break;
      }
    } else {
      breakIterationCounter = 0;
    }
  }

  xOpt.resize(d);
  xOpt = x;
  fOpt = fx;
  Printer::getInstance().printStatusEnd();
}

ScalarFunctionGradient& AdaptiveGradientDescent::getObjectiveGradient() const { return *fGradient; }

double AdaptiveGradientDescent::getTolerance() const { return theta; }

void AdaptiveGradientDescent::setTolerance(double tolerance) { theta = tolerance; }

double AdaptiveGradientDescent::getStepSizeIncreaseFactor() const { return rhoAlphaPlus; }

void AdaptiveGradientDescent::setStepSizeIncreaseFactor(double stepSizeIncreaseFactor) {
  rhoAlphaPlus = stepSizeIncreaseFactor;
}

double AdaptiveGradientDescent::getStepSizeDecreaseFactor() const { return rhoAlphaMinus; }

void AdaptiveGradientDescent::setStepSizeDecreaseFactor(double stepSizeDecreaseFactor) {
  rhoAlphaMinus = stepSizeDecreaseFactor;
}

double AdaptiveGradientDescent::getLineSearchAccuracy() const { return rhoLs; }

void AdaptiveGradientDescent::setLineSearchAccuracy(double lineSearchAccuracy) {
  rhoLs = lineSearchAccuracy;
}

void AdaptiveGradientDescent::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new AdaptiveGradientDescent(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

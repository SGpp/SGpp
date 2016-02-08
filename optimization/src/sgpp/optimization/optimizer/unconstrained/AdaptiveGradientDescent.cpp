// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>

namespace SGPP {
namespace optimization {
namespace optimizer {

AdaptiveGradientDescent::AdaptiveGradientDescent(
  ScalarFunction& f,
  ScalarFunctionGradient& fGradient,
  size_t maxItCount,
  float_t tolerance,
  float_t stepSizeIncreaseFactor,
  float_t stepSizeDecreaseFactor,
  float_t lineSearchAccuracy) :
  UnconstrainedOptimizer(f, maxItCount),
  fGradient(fGradient),
  theta(tolerance),
  rhoAlphaPlus(stepSizeIncreaseFactor),
  rhoAlphaMinus(stepSizeDecreaseFactor),
  rhoLs(lineSearchAccuracy) {
}

AdaptiveGradientDescent::~AdaptiveGradientDescent() {
}

void AdaptiveGradientDescent::optimize() {
  Printer::getInstance().printStatusBegin("Optimizing (adaptive gradient descent)...");

  const size_t d = f.getNumberOfParameters();

  xOpt.resize(0);
  fOpt = NAN;
  xHist.resize(0, d);
  fHist.resize(0);

  base::DataVector x(x0);
  float_t fx = f.eval(x);
  base::DataVector gradFx(d);

  xHist.appendRow(x);
  fHist.append(fx);

  base::DataVector xNew(d);
  float_t fxNew;

  size_t k = 0;
  float_t alpha = 1.0;
  base::DataVector dir(d);
  bool inDomain;

  size_t breakIterationCounter = 0;
  const size_t BREAK_ITERATION_COUNTER_MAX = 10;

  while (k < N) {
    // calculate gradient and norm
    fx = fGradient.eval(x, gradFx);
    k++;

    const float_t gradFxNorm = gradFx.l2Norm();

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
    fxNew = (inDomain ? f.eval(xNew) : INFINITY);
    k++;

    // inner product of gradient and search direction
    const float_t gradFxTimesDir = -gradFxNorm;

    // line search
    while ((fxNew > fx + rhoLs * alpha * gradFxTimesDir) &&
           (alpha > 0.0)) {
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
      fxNew = (inDomain ? f.eval(xNew) : INFINITY);
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
    Printer::getInstance().printStatusUpdate(
      std::to_string(k) + " evaluations, x = " + x.toString() +
      ", f(x) = " + std::to_string(fx));

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

ScalarFunctionGradient& AdaptiveGradientDescent::getObjectiveGradient() const {
  return fGradient;
}

float_t AdaptiveGradientDescent::getTolerance() const {
  return theta;
}

void AdaptiveGradientDescent::setTolerance(float_t tolerance) {
  theta = tolerance;
}

float_t AdaptiveGradientDescent::getStepSizeIncreaseFactor() const {
  return rhoAlphaPlus;
}

void AdaptiveGradientDescent::setStepSizeIncreaseFactor(
  float_t stepSizeIncreaseFactor) {
  rhoAlphaPlus = stepSizeIncreaseFactor;
}

float_t AdaptiveGradientDescent::getStepSizeDecreaseFactor() const {
  return rhoAlphaMinus;
}

void AdaptiveGradientDescent::setStepSizeDecreaseFactor(
  float_t stepSizeDecreaseFactor) {
  rhoAlphaMinus = stepSizeDecreaseFactor;
}

float_t AdaptiveGradientDescent::getLineSearchAccuracy() const {
  return rhoLs;
}

void AdaptiveGradientDescent::setLineSearchAccuracy(
  float_t lineSearchAccuracy) {
  rhoLs = lineSearchAccuracy;
}

void AdaptiveGradientDescent::clone(
  std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(
            new AdaptiveGradientDescent(*this));
}
}
}
}

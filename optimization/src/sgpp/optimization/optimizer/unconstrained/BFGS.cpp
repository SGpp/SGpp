// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/BFGS.hpp>

#include <limits>

namespace sgpp {
namespace optimization {
namespace optimizer {

BFGS::BFGS(const base::ScalarFunction& f, const base::ScalarFunctionGradient& fGradient,
           size_t maxItCount, double tolerance, double stepSizeIncreaseFactor,
           double stepSizeDecreaseFactor, double lineSearchAccuracy)
    : UnconstrainedOptimizer(f, &fGradient, nullptr, maxItCount),
      theta(tolerance),
      rhoAlphaPlus(stepSizeIncreaseFactor),
      rhoAlphaMinus(stepSizeDecreaseFactor),
      rhoLs(lineSearchAccuracy) {
}

BFGS::BFGS(const BFGS& other)
    : UnconstrainedOptimizer(other),
      theta(other.theta),
      rhoAlphaPlus(other.rhoAlphaPlus),
      rhoAlphaMinus(other.rhoAlphaMinus),
      rhoLs(other.rhoLs) {
}

BFGS::~BFGS() {}

void BFGS::optimize() {
  base::Printer::getInstance().printStatusBegin("Optimizing (BFGS)...");

  const size_t d = f->getNumberOfParameters();

  xOpt.resize(0);
  fOpt = std::numeric_limits<double>::quiet_NaN();
  xHist.resize(0, d);
  fHist.resize(0);

  base::DataVector x(x0);
  double fx = std::numeric_limits<double>::quiet_NaN();
  base::DataVector gradFx(d);

  base::DataVector xNew(d);
  double fxNew;
  base::DataVector gradFxNew(d);

  base::DataVector delta(d);
  base::DataVector y(d);

  base::DataMatrix inverseHessian(d, d);
  base::DataMatrix inverseHessianNew(d, d);
  base::DataMatrix M(d, d);

  for (size_t i = 0; i < d; i++) {
    for (size_t j = 0; j < d; j++) {
      inverseHessian(i, j) = (i == j ? 1.0 : 0.0);
    }
  }

  size_t k = 0;
  double alpha = 1.0;
  base::DataVector dir(d);
  bool inDomain;

  size_t breakIterationCounter = 0;
  const size_t BREAK_ITERATION_COUNTER_MAX = 10;

  while (k < N) {
    // calculate gradient
    fx = fGradient->eval(x, gradFx);
    k++;

    if (k == 1) {
      xHist.appendRow(x);
      fHist.append(fx);
    }

    const double gradFxNorm = gradFx.l2Norm();

    if (gradFxNorm == 0.0) {
      break;
    }

    for (size_t i = 0; i < d; i++) {
      dir[i] = 0.0;

      for (size_t j = 0; j < d; j++) {
        dir[i] -= inverseHessian(i, j) * gradFx[j];
      }
    }

    if (dir.dotProduct(gradFx) > 0.0) {
      for (size_t t = 0; t < d; t++) {
        dir[t] = -gradFx[t] / gradFxNorm;
      }
    }

    inDomain = true;

    for (size_t t = 0; t < d; t++) {
      xNew[t] = x[t] + alpha * dir[t];

      if ((xNew[t] < 0.0) || (xNew[t] > 1.0)) {
        inDomain = false;
        break;
      }
    }

    // evaluate at new point
    fxNew = (inDomain ? f->eval(xNew) : std::numeric_limits<double>::infinity());
    k++;

    // inner product of gradient and search direction
    const double gradFxTimesDir = gradFx.dotProduct(dir);

    // line search
    while (fxNew > fx + rhoLs * alpha * gradFxTimesDir) {
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
      fxNew = (inDomain ? fGradient->eval(xNew, gradFxNew) :
          std::numeric_limits<double>::infinity());
      k++;
    }

    for (size_t t = 0; t < d; t++) {
      delta[t] = alpha * dir[t];
      y[t] = gradFxNew[t] - gradFx[t];
    }

    // save new point
    x = xNew;
    fx = fxNew;
    gradFx = gradFxNew;
    xHist.appendRow(x);
    fHist.append(fx);

    // increase step size
    alpha *= rhoAlphaPlus;

    const double deltaTimesY = delta.dotProduct(y);

    if (deltaTimesY != 0.0) {
      for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
          M(i, j) = (i == j ? 1.0 : 0.0) - y[i] * delta[j] / deltaTimesY;
        }
      }

      for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
          double entry = delta[i] * delta[j] / deltaTimesY;

          for (size_t p = 0; p < d; p++) {
            for (size_t q = 0; q < d; q++) {
              entry += M(p, i) * inverseHessian(p, q) * M(q, j);
            }
          }

          inverseHessianNew(i, j) = entry;
        }
      }

      inverseHessian = inverseHessianNew;
    }

    // status printing
    base::Printer::getInstance().printStatusUpdate(
        std::to_string(k) + " evaluations, x = " + x.toString() + ", f(x) = " + std::to_string(fx));

    // stopping criterion:
    // stop if delta is smaller than tolerance theta
    // in BREAK_ITERATION_COUNTER_MAX consecutive iterations
    if (delta.l2Norm() < theta) {
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
  base::Printer::getInstance().printStatusEnd();
}

double BFGS::getTolerance() const { return theta; }

void BFGS::setTolerance(double tolerance) { theta = tolerance; }

double BFGS::getStepSizeIncreaseFactor() const { return rhoAlphaPlus; }

void BFGS::setStepSizeIncreaseFactor(double stepSizeIncreaseFactor) {
  rhoAlphaPlus = stepSizeIncreaseFactor;
}

double BFGS::getStepSizeDecreaseFactor() const { return rhoAlphaMinus; }

void BFGS::setStepSizeDecreaseFactor(double stepSizeDecreaseFactor) {
  rhoAlphaMinus = stepSizeDecreaseFactor;
}

double BFGS::getLineSearchAccuracy() const { return rhoLs; }

void BFGS::setLineSearchAccuracy(double lineSearchAccuracy) { rhoLs = lineSearchAccuracy; }

void BFGS::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new BFGS(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

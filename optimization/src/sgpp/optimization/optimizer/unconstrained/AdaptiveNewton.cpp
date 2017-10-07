// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveNewton.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>

#include <algorithm>

namespace sgpp {
namespace optimization {
namespace optimizer {

AdaptiveNewton::AdaptiveNewton(const ScalarFunction& f, const ScalarFunctionHessian& fHessian,
                               size_t maxItCount, double tolerance, double stepSizeIncreaseFactor,
                               double stepSizeDecreaseFactor, double dampingIncreaseFactor,
                               double dampingDecreaseFactor, double lineSearchAccuracy)
    : UnconstrainedOptimizer(f, nullptr, &fHessian, maxItCount),
      theta(tolerance),
      rhoAlphaPlus(stepSizeIncreaseFactor),
      rhoAlphaMinus(stepSizeDecreaseFactor),
      rhoLambdaPlus(dampingIncreaseFactor),
      rhoLambdaMinus(dampingDecreaseFactor),
      rhoLs(lineSearchAccuracy),
      defaultSleSolver(sle_solver::GaussianElimination()),
      sleSolver(defaultSleSolver) {
}

AdaptiveNewton::AdaptiveNewton(const ScalarFunction& f, const ScalarFunctionHessian& fHessian,
                               size_t maxItCount, double tolerance, double stepSizeIncreaseFactor,
                               double stepSizeDecreaseFactor, double dampingIncreaseFactor,
                               double dampingDecreaseFactor, double lineSearchAccuracy,
                               const sle_solver::SLESolver& sleSolver)
    : UnconstrainedOptimizer(f, nullptr, &fHessian, maxItCount),
      theta(tolerance),
      rhoAlphaPlus(stepSizeIncreaseFactor),
      rhoAlphaMinus(stepSizeDecreaseFactor),
      rhoLambdaPlus(dampingIncreaseFactor),
      rhoLambdaMinus(dampingDecreaseFactor),
      rhoLs(lineSearchAccuracy),
      defaultSleSolver(sle_solver::GaussianElimination()),
      sleSolver(sleSolver) {
}

AdaptiveNewton::AdaptiveNewton(const AdaptiveNewton& other)
    : UnconstrainedOptimizer(other),
      theta(other.theta),
      rhoAlphaPlus(other.rhoAlphaPlus),
      rhoAlphaMinus(other.rhoAlphaMinus),
      rhoLambdaPlus(other.rhoLambdaPlus),
      rhoLambdaMinus(other.rhoLambdaMinus),
      rhoLs(other.rhoLs),
      defaultSleSolver(sle_solver::GaussianElimination()),
      sleSolver(other.sleSolver) {
}

AdaptiveNewton::~AdaptiveNewton() {}

void AdaptiveNewton::optimize() {
  Printer::getInstance().printStatusBegin("Optimizing (adaptive Newton)...");

  const size_t d = f->getNumberOfParameters();

  xOpt.resize(0);
  fOpt = NAN;
  xHist.resize(0, d);
  fHist.resize(0);

  base::DataVector x(x0);
  double fx = NAN;
  base::DataVector gradFx(d);
  base::DataMatrix hessianFx(d, d);

  base::DataVector b(d);
  bool lsSolved;

  base::DataVector xNew(x0);
  double fxNew;

  FullSLE system(hessianFx);
  size_t k = 0;
  double alpha = 1.0;
  double lambda = 1.0;
  base::DataVector dir(d);
  bool inDomain;

  size_t breakIterationCounter = 0;
  const size_t BREAK_ITERATION_COUNTER_MAX = 10;

  const double ALPHA1 = 1e-6;
  const double ALPHA2 = 1e-6;
  const double P = 0.1;
  const bool statusPrintingEnabled = Printer::getInstance().isStatusPrintingEnabled();

  while (k < N) {
    // calculate gradient and Hessian
    fx = fHessian->eval(x, gradFx, hessianFx);
    k++;

    if (k == 1) {
      xHist.appendRow(x);
      fHist.append(fx);
    }

    const double gradFxNorm = gradFx.l2Norm();

    if (gradFxNorm == 0.0) {
      break;
    }

    for (size_t t = 0; t < d; t++) {
      // RHS of linear system to be solved
      b[t] = -gradFx[t];
      // add damping
      hessianFx(t, t) += lambda;
    }

    // solve linear system with damped Hessian as system matrix
    if (statusPrintingEnabled) {
      Printer::getInstance().disableStatusPrinting();
    }

    lsSolved = sleSolver.solve(system, b, dir);

    if (statusPrintingEnabled) {
      Printer::getInstance().enableStatusPrinting();
    }

    const double dirNorm = dir.l2Norm();

    // acceptance criterion
    if (lsSolved && (b.dotProduct(dir) >=
                     std::min(ALPHA1, ALPHA2 * std::pow(dirNorm, P)) * dirNorm * dirNorm)) {
      // normalize search direction
      for (size_t t = 0; t < d; t++) {
        dir[t] /= dirNorm;
      }
    } else {
      // restart method
      // (negated normalized gradient as new search direction)
      for (size_t t = 0; t < d; t++) {
        dir[t] = b[t] / gradFxNorm;
      }
    }

    inDomain = true;

    for (size_t t = 0; t < d; t++) {
      // new point
      xNew[t] = x[t] + alpha * dir[t];

      if ((xNew[t] < 0.0) || (xNew[t] > 1.0)) {
        inDomain = false;
        break;
      }
    }

    // evaluate at new point
    fxNew = (inDomain ? f->eval(xNew) : INFINITY);
    k++;

    // inner product of gradient and search direction
    double gradFxTimesDir = gradFx.dotProduct(dir);

    // line search
    while ((fxNew > fx + rhoLs * alpha * gradFxTimesDir) && (alpha > 0.0)) {
      alpha *= rhoAlphaMinus;

      // increase damping
      if (rhoLambdaPlus != 1.0) {
        const double oldLambda = lambda;
        lambda *= rhoLambdaPlus;

        for (size_t t = 0; t < d; t++) {
          // add damping
          hessianFx(t, t) += lambda - oldLambda;
        }

        // solve linear system with damped Hessian as system matrix
        Printer::getInstance().disableStatusPrinting();
        lsSolved = sleSolver.solve(system, b, dir);
        Printer::getInstance().enableStatusPrinting();

        // recalculate inner product
        gradFxTimesDir = gradFx.dotProduct(dir);
      }

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
    alpha = std::min(rhoAlphaPlus * alpha, 1.0);

    // decrease damping
    lambda *= rhoLambdaMinus;

    // status printing
    Printer::getInstance().printStatusUpdate(std::to_string(k) + " evaluations, x = " +
                                             x.toString() + ", f(x) = " + std::to_string(fx));

    // stopping criterion:
    // stop if alpha * dir is smaller than tolerance theta
    // in BREAK_ITERATION_COUNTER_MAX consecutive iterations
    if (alpha * dir.l2Norm() < theta) {
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

double AdaptiveNewton::getTolerance() const { return theta; }

void AdaptiveNewton::setTolerance(double tolerance) { theta = tolerance; }

double AdaptiveNewton::getStepSizeIncreaseFactor() const { return rhoAlphaPlus; }

void AdaptiveNewton::setStepSizeIncreaseFactor(double stepSizeIncreaseFactor) {
  rhoAlphaPlus = stepSizeIncreaseFactor;
}

double AdaptiveNewton::getStepSizeDecreaseFactor() const { return rhoAlphaMinus; }

void AdaptiveNewton::setStepSizeDecreaseFactor(double stepSizeDecreaseFactor) {
  rhoAlphaMinus = stepSizeDecreaseFactor;
}

double AdaptiveNewton::getDampingIncreaseFactor() const { return rhoLambdaPlus; }

void AdaptiveNewton::setDampingIncreaseFactor(double dampingIncreaseFactor) {
  rhoLambdaPlus = dampingIncreaseFactor;
}

double AdaptiveNewton::getDampingDecreaseFactor() const { return rhoLambdaMinus; }

void AdaptiveNewton::setDampingDecreaseFactor(double dampingDecreaseFactor) {
  rhoLambdaMinus = dampingDecreaseFactor;
}

double AdaptiveNewton::getLineSearchAccuracy() const { return rhoLs; }

void AdaptiveNewton::setLineSearchAccuracy(double lineSearchAccuracy) {
  rhoLs = lineSearchAccuracy;
}

void AdaptiveNewton::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new AdaptiveNewton(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

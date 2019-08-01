// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/system/FullSLE.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveNewton.hpp>

#include <algorithm>
#include <vector>

namespace sgpp {
namespace optimization {
namespace optimizer {

AdaptiveNewton::AdaptiveNewton(const base::ScalarFunction& f,
                               const base::ScalarFunctionHessian& fHessian, size_t maxItCount,
                               double tolerance, double stepSizeIncreaseFactor,
                               double stepSizeDecreaseFactor, double dampingIncreaseFactor,
                               double dampingDecreaseFactor, double lineSearchAccuracy)
    : UnconstrainedOptimizer(f, nullptr, &fHessian, maxItCount),
      theta(tolerance),
      rhoAlphaPlus(stepSizeIncreaseFactor),
      rhoAlphaMinus(stepSizeDecreaseFactor),
      rhoLambdaPlus(dampingIncreaseFactor),
      rhoLambdaMinus(dampingDecreaseFactor),
      rhoLs(lineSearchAccuracy),
      defaultSleSolver(base::sle_solver::GaussianElimination()),
      sleSolver(defaultSleSolver) {
}

AdaptiveNewton::AdaptiveNewton(const base::ScalarFunction& f,
                               const base::ScalarFunctionHessian& fHessian, size_t maxItCount,
                               double tolerance, double stepSizeIncreaseFactor,
                               double stepSizeDecreaseFactor, double dampingIncreaseFactor,
                               double dampingDecreaseFactor, double lineSearchAccuracy,
                               const base::sle_solver::SLESolver& sleSolver)
    : UnconstrainedOptimizer(f, maxItCount),
      theta(tolerance),
      rhoAlphaPlus(stepSizeIncreaseFactor),
      rhoAlphaMinus(stepSizeDecreaseFactor),
      rhoLambdaPlus(dampingIncreaseFactor),
      rhoLambdaMinus(dampingDecreaseFactor),
      rhoLs(lineSearchAccuracy),
      defaultSleSolver(base::sle_solver::GaussianElimination()),
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
      defaultSleSolver(base::sle_solver::GaussianElimination()),
      sleSolver(other.sleSolver) {
}

AdaptiveNewton::~AdaptiveNewton() {}

void AdaptiveNewton::optimize() {
  base::Printer::getInstance().printStatusBegin("Optimizing (adaptive Newton)...");

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

  base::FullSLE system(hessianFx);
  size_t k = 0;
  double alpha = 1.0;
  double lambda = 1.0;
  base::DataVector dir(d);
  bool inDomain;
  std::vector<bool> isConstraintBinding(d, false);

  size_t breakIterationCounter = 0;
  const size_t BREAK_ITERATION_COUNTER_MAX = 10;

  const double ALPHA1 = 1e-6;
  const double ALPHA2 = 1e-6;
  const double BINDING_TOLERANCE = 1e-6;
  const double P = 0.1;
  const bool statusPrintingEnabled = base::Printer::getInstance().isStatusPrintingEnabled();

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

      // eliminate rows and columns corresponding to binding constraints
      if (isConstraintBinding[t]) {
        // discard variable by setting RHS to zero
        b[t] = 0.0;
        // eliminate variable from linear system matrix
        hessianFx(t, t) = 1.0;

        for (size_t t2 = 0; t2 < d; t2++) {
          if (t != t2) {
            hessianFx(t, t2) = 0.0;
            hessianFx(t2, t) = 0.0;
          }
        }
      }
    }

    // solve linear system with damped Hessian as system matrix
    if (statusPrintingEnabled) {
      base::Printer::getInstance().disableStatusPrinting();
    }

    lsSolved = sleSolver.solve(system, b, dir);

    if (statusPrintingEnabled) {
      base::Printer::getInstance().enableStatusPrinting();
    }

    double dirNorm = dir.l2Norm();

    // acceptance criterion
    if (!(lsSolved && (b.dotProduct(dir) >=
                       std::min(ALPHA1, ALPHA2 * std::pow(dirNorm, P)) * dirNorm * dirNorm))) {
      // restart method (negated normalized gradient as new search direction)
      dir = b;
    }

      for (size_t t = 0; t < d; t++) {
      // is constraint binding?
      // (i.e., we are at the boundary and the search direction points outwards)
      if (((x[t] < BINDING_TOLERANCE) && (dir[t] < 0.0)) ||
          ((x[t] > 1.0 - BINDING_TOLERANCE) && (dir[t] > 0.0))) {
        // discard variable by setting direction to zero
        dir[t] = 0.0;

        // was the constraint not binding in the previous iteration?
        if (!isConstraintBinding[t]) {
          // reset step size as it's most likely very small due to approach to the boundary
          alpha = 1.0;
          isConstraintBinding[t] = true;
      }
    } else {
        isConstraintBinding[t] = false;
      }
    }

    dirNorm = dir.l2Norm();

    if (dirNorm == 0.0) {
      break;
    }

    inDomain = true;

    for (size_t t = 0; t < d; t++) {
      // normalize search direction
      dir[t] /= dirNorm;

      // new point
      xNew[t] = x[t] + alpha * dir[t];

      if ((xNew[t] < 0.0) || (xNew[t] > 1.0)) {
        inDomain = false;
        break;
      }
    }

    // evaluate at new point
    if (inDomain) {
      fxNew = f->eval(xNew);
    k++;
    } else {
      fxNew = INFINITY;
    }

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

          // eliminate rows and columns corresponding to binding constraints
          if (isConstraintBinding[t]) {
            // discard variable by setting RHS to zero
            b[t] = 0.0;
            // eliminate variable from linear system matrix
            hessianFx(t, t) = 1.0;

            for (size_t t2 = 0; t2 < d; t2++) {
              if (t != t2) {
                hessianFx(t, t2) = 0.0;
                hessianFx(t2, t) = 0.0;
              }
            }
          }
        }

        // solve linear system with damped Hessian as system matrix
        if (statusPrintingEnabled) {
          base::Printer::getInstance().disableStatusPrinting();
        }

        lsSolved = sleSolver.solve(system, b, dir);

        if (statusPrintingEnabled) {
          base::Printer::getInstance().enableStatusPrinting();
        }

        for (size_t t = 0; t < d; t++) {
          // is constraint binding?
          // (i.e., we are at the boundary and the search direction points outwards)
          if (((x[t] < BINDING_TOLERANCE) && (dir[t] < 0.0)) ||
              ((x[t] > 1.0 - BINDING_TOLERANCE) && (dir[t] > 0.0))) {
            // discard variable by setting direction to zero
            dir[t] = 0.0;

            // was the constraint not binding in the previous iteration?
            if (!isConstraintBinding[t]) {
              // reset step size as it's most likely very small due to approach to the boundary
              alpha = 1.0;
              isConstraintBinding[t] = true;
            }
          } else {
            isConstraintBinding[t] = false;
          }
        }

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
      if (inDomain) {
        fxNew = f->eval(xNew);
        k++;
      } else {
        fxNew = INFINITY;
      }
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
    base::Printer::getInstance().printStatusUpdate(
        std::to_string(k) + " evaluations, x = " + x.toString() + ", f(x) = " + std::to_string(fx));

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
  base::Printer::getInstance().printStatusEnd();
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

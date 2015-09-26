// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveNewton.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      AdaptiveNewton::AdaptiveNewton(
        ScalarFunction& f,
        ScalarFunctionHessian& fHessian,
        size_t maxItCount,
        float_t tolerance,
        float_t stepSizeIncreaseFactor,
        float_t stepSizeDecreaseFactor,
        float_t dampingIncreaseFactor,
        float_t dampingDecreaseFactor,
        float_t lineSearchAccuracy) :
        UnconstrainedOptimizer(f, maxItCount),
        fHessian(fHessian),
        theta(tolerance),
        rhoAlphaPlus(stepSizeIncreaseFactor),
        rhoAlphaMinus(stepSizeDecreaseFactor),
        rhoLambdaPlus(dampingIncreaseFactor),
        rhoLambdaMinus(dampingDecreaseFactor),
        rhoLs(lineSearchAccuracy),
        defaultSleSolver(sle_solver::GaussianElimination()),
        sleSolver(defaultSleSolver) {
      }

      AdaptiveNewton::AdaptiveNewton(
        ScalarFunction& f,
        ScalarFunctionHessian& fHessian,
        size_t maxItCount,
        float_t tolerance,
        float_t stepSizeIncreaseFactor,
        float_t stepSizeDecreaseFactor,
        float_t dampingIncreaseFactor,
        float_t dampingDecreaseFactor,
        float_t lineSearchAccuracy,
        const sle_solver::SLESolver& sleSolver) :
        UnconstrainedOptimizer(f, N),
        fHessian(fHessian),
        theta(tolerance),
        rhoAlphaPlus(stepSizeIncreaseFactor),
        rhoAlphaMinus(stepSizeDecreaseFactor),
        rhoLambdaPlus(dampingIncreaseFactor),
        rhoLambdaMinus(dampingDecreaseFactor),
        rhoLs(lineSearchAccuracy),
        defaultSleSolver(sle_solver::GaussianElimination()),
        sleSolver(sleSolver) {
      }

      float_t AdaptiveNewton::optimize(base::DataVector& xOpt) {
        printer.printStatusBegin("Optimizing (adaptive Newton)...");

        const size_t d = f.getDimension();
        base::DataVector x(x0);
        float_t fx = NAN;
        base::DataVector gradFx(d);
        base::DataMatrix hessianFx(d, d);

        base::DataVector b(d);
        bool lsSolved;

        base::DataVector xNew(x0);
        float_t fxNew;

        FullSLE system(hessianFx);
        size_t k = 0;
        float_t alpha = 1.0;
        float_t lambda = 1.0;
        base::DataVector dir(d);
        bool inDomain;

        size_t breakIterationCounter = 0;
        const size_t BREAK_ITERATION_COUNTER_MAX = 10;

        const float_t ALPHA1 = 1e-6;
        const float_t ALPHA2 = 1e-6;
        const float_t P = 0.1;
        const bool statusPrintingEnabled = printer.isStatusPrintingEnabled();

        while (k < N) {
          // calculate gradient and Hessian
          fx = fHessian.eval(x, gradFx, hessianFx);
          k++;

          const float_t gradFxNorm = gradFx.l2Norm();

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
          system.setA(hessianFx);

          if (statusPrintingEnabled) {
            printer.disableStatusPrinting();
          }

          lsSolved = sleSolver.solve(system, b, dir);

          if (statusPrintingEnabled) {
            printer.enableStatusPrinting();
          }

          const float_t dirNorm = dir.l2Norm();

          // acceptance criterion
          if (lsSolved && (b.dotProduct(dir) >=
                           std::min(ALPHA1, ALPHA2 * std::pow(dirNorm, P)) *
                           dirNorm * dirNorm)) {
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
          fxNew = (inDomain ? f.eval(xNew) : INFINITY);
          k++;

          // inner product of gradient and search direction
          float_t gradFxTimesDir = gradFx.dotProduct(dir);

          // line search
          while (fxNew > fx + rhoLs * alpha * gradFxTimesDir) {
            alpha *= rhoAlphaMinus;

            // increase damping
            if (rhoLambdaPlus != 1.0) {
              const float_t oldLambda = lambda;
              lambda *= rhoLambdaPlus;

              for (size_t t = 0; t < d; t++) {
                // add damping
                hessianFx(t, t) += lambda - oldLambda;
              }

              // solve linear system with damped Hessian as system matrix
              system.setA(hessianFx);
              printer.disableStatusPrinting();
              lsSolved = sleSolver.solve(system, b, dir);
              printer.enableStatusPrinting();

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
            fxNew = (inDomain ? f.eval(xNew) : INFINITY);
            k++;
          }

          // save new point
          x = xNew;
          fx = fxNew;

          // increase step size
          alpha = std::min(rhoAlphaPlus * alpha, float_t(1.0) );

          // decrease damping
          lambda *= rhoLambdaMinus;

          // status printing
          printer.printStatusUpdate(
            std::to_string(k) + " evaluations, x = " + x.toString() +
            ", f(x) = " + std::to_string(fx));

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
        printer.printStatusEnd();

        return fx;
      }

      ScalarFunctionHessian& AdaptiveNewton::getObjectiveHessian() const {
        return fHessian;
      }

      float_t AdaptiveNewton::getTolerance() const {
        return theta;
      }

      void AdaptiveNewton::setTolerance(float_t tolerance) {
        theta = tolerance;
      }

      float_t AdaptiveNewton::getStepSizeIncreaseFactor() const {
        return rhoAlphaPlus;
      }

      void AdaptiveNewton::setStepSizeIncreaseFactor(
        float_t stepSizeIncreaseFactor) {
        rhoAlphaPlus = stepSizeIncreaseFactor;
      }

      float_t AdaptiveNewton::getStepSizeDecreaseFactor() const {
        return rhoAlphaMinus;
      }

      void AdaptiveNewton::setStepSizeDecreaseFactor(
        float_t stepSizeDecreaseFactor) {
        rhoAlphaMinus = stepSizeDecreaseFactor;
      }

      float_t AdaptiveNewton::getDampingIncreaseFactor() const {
        return rhoLambdaPlus;
      }

      void AdaptiveNewton::setDampingIncreaseFactor(
        float_t dampingIncreaseFactor) {
        rhoLambdaPlus = dampingIncreaseFactor;
      }

      float_t AdaptiveNewton::getDampingDecreaseFactor() const {
        return rhoLambdaMinus;
      }

      void AdaptiveNewton::setDampingDecreaseFactor(
        float_t dampingDecreaseFactor) {
        rhoLambdaMinus = dampingDecreaseFactor;
      }

      float_t AdaptiveNewton::getLineSearchAccuracy() const {
        return rhoLs;
      }

      void AdaptiveNewton::setLineSearchAccuracy(
        float_t lineSearchAccuracy) {
        rhoLs = lineSearchAccuracy;
      }

    }
  }
}

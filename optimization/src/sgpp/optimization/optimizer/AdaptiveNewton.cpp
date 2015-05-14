// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/AdaptiveNewton.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      AdaptiveNewton::AdaptiveNewton(
        ObjectiveFunction& f,
        ObjectiveHessian& fHessian,
        size_t maxItCount,
        float_t tolerance,
        float_t stepSizeIncreaseFactor,
        float_t stepSizeDecreaseFactor,
        float_t dampingIncreaseFactor,
        float_t dampingDecreaseFactor,
        float_t lineSearchAccuracy) :
        Optimizer(f, N),
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
        ObjectiveFunction& f,
        ObjectiveHessian& fHessian,
        size_t maxItCount,
        float_t tolerance,
        float_t stepSizeIncreaseFactor,
        float_t stepSizeDecreaseFactor,
        float_t dampingIncreaseFactor,
        float_t dampingDecreaseFactor,
        float_t lineSearchAccuracy,
        const sle_solver::SLESolver& sleSolver) :
        Optimizer(f, N),
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
        size_t k;
        float_t alpha = 1.0;
        float_t lambda = 1.0;
        base::DataVector dir(d);

        size_t breakIterationCounter = 0;
        const size_t BREAK_ITERATION_COUNTER_MAX = 10;

        for (k = 0; k < N; k++) {
          // calculate gradient, Hessian and gradient norm
          fx = fHessian.eval(x, gradFx, hessianFx);

          // DEBUG
          /*std::cout << "\nk = " << k << "\n";
          std::cout << "x = " << x.toString() << "\n";
          std::cout << "fx = " << fx << "\n";
          std::cout << "gradFx = " << gradFx.toString() << "\n";
          std::cout << "hessianFx = " << hessianFx.toString() << "\n";*/

          for (size_t t = 0; t < d; t++) {
            // RHS of linear system to be solved
            b[t] = -gradFx[t];
            // add damping
            hessianFx.set(t, t, hessianFx.get(t, t) + lambda);
          }

          // solve linear system with damped Hessian as system matrix
          system.setA(hessianFx);
          printer.disableStatusPrinting();
          lsSolved = sleSolver.solve(system, b, dir);
          printer.enableStatusPrinting();

          if (!lsSolved) {
            // restart method
            // (negated normalized gradient as new search direction)
            const float_t gradFxNorm = gradFx.l2Norm();

            for (size_t t = 0; t < d; t++) {
              dir[t] = gradFx[t] / gradFxNorm;
            }
          }

          for (size_t t = 0; t < d; t++) {
            // new point
            xNew[t] = x[t] + alpha * dir[t];
          }

          // evaluate at new point
          fxNew = f.eval(xNew);

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
                hessianFx.set(t, t, hessianFx.get(t, t) - oldLambda + lambda);
              }

              // solve linear system with damped Hessian as system matrix
              system.setA(hessianFx);
              printer.disableStatusPrinting();
              lsSolved = sleSolver.solve(system, b, dir);
              printer.enableStatusPrinting();

              // recalculate inner product
              gradFxTimesDir = gradFx.dotProduct(dir);
            }

            // recalculate new point
            for (size_t t = 0; t < d; t++) {
              xNew[t] = x[t] + alpha * dir[t];
            }

            // evaluate at new point
            fxNew = f.eval(xNew);
          }

          // save new point
          x = xNew;
          fx = fxNew;

          // increase step size
          alpha = std::min(rhoAlphaPlus * alpha, 1.0);

          // decrease damping
          lambda *= rhoLambdaMinus;

          // status printing
          printer.printStatusUpdate(std::to_string(k) + " steps, f(x) = " +
                                    std::to_string(fx));

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

        printer.printStatusUpdate(std::to_string(k) + " steps, f(x) = " +
                                  std::to_string(fx));
        printer.printStatusEnd();

        return fx;
      }

      ObjectiveHessian& AdaptiveNewton::getObjectiveHessian() const {
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

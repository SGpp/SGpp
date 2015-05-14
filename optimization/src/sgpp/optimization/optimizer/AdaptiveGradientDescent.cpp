// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/AdaptiveGradientDescent.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      AdaptiveGradientDescent::AdaptiveGradientDescent(
        ObjectiveFunction& f,
        ObjectiveGradient& fGradient,
        size_t maxItCount,
        float_t tolerance,
        float_t stepSizeIncreaseFactor,
        float_t stepSizeDecreaseFactor,
        float_t lineSearchAccuracy) :
        Optimizer(f, N),
        fGradient(fGradient),
        theta(tolerance),
        rhoAlphaPlus(stepSizeIncreaseFactor),
        rhoAlphaMinus(stepSizeDecreaseFactor),
        rhoLs(lineSearchAccuracy) {
      }

      float_t AdaptiveGradientDescent::optimize(base::DataVector& xOpt) {
        printer.printStatusBegin("Optimizing (adaptive gradient descent)...");

        const size_t d = f.getDimension();

        base::DataVector x(x0);
        float_t fx;
        base::DataVector gradFx(d);

        base::DataVector xNew(d);
        float_t fxNew;

        size_t k;
        float_t alpha = 1.0;
        base::DataVector dir(d);

        size_t breakIterationCounter = 0;
        const size_t BREAK_ITERATION_COUNTER_MAX = 10;

        for (k = 0; k < N; k++) {
          // calculate gradient and norm
          fx = fGradient.eval(x, gradFx);
          const float_t gradFxNorm = gradFx.l2Norm();

          for (size_t t = 0; t < d; t++) {
            // search direction (normalized negated gradient)
            dir[t] = -gradFx[t] / gradFxNorm;
            // new point
            xNew[t] = x[t] + alpha * dir[t];
          }

          // evaluate at new point
          fxNew = f.eval(xNew);

          // inner product of gradient and search direction
          const float_t gradFxTimesDir = -gradFxNorm;

          // line search
          while (fxNew > fx + rhoLs * alpha * gradFxTimesDir) {
            alpha *= rhoAlphaMinus;

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
          alpha *= rhoAlphaPlus;

          // status printing
          printer.printStatusUpdate(std::to_string(k) + " steps, f(x) = " +
                                    std::to_string(fx));

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

        printer.printStatusUpdate(std::to_string(k) + " steps, f(x) = " +
                                  std::to_string(fx));
        printer.printStatusEnd();

        return fx;
      }

      ObjectiveGradient& AdaptiveGradientDescent::getObjectiveGradient() const {
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

    }
  }
}

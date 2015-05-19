// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/BFGS.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      BFGS::BFGS(
        ObjectiveFunction& f,
        ObjectiveGradient& fGradient,
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

      float_t BFGS::optimize(base::DataVector& xOpt) {
        printer.printStatusBegin("Optimizing (BFGS)...");

        const size_t d = f.getDimension();

        base::DataVector x(x0);
        float_t fx = NAN;
        base::DataVector gradFx(d);

        base::DataVector xNew(d);
        float_t fxNew;
        base::DataVector gradFxNew(d);

        base::DataVector delta(d);
        base::DataVector y(d);

        base::DataMatrix inverseHessian(d, d);
        base::DataMatrix inverseHessianNew(d, d);
        base::DataMatrix M(d, d);

        for (size_t i = 0; i < d; i++) {
          for (size_t j = 0; j < d; j++) {
            inverseHessian.set(i, j, (i == j ? 1.0 : 0.0));
          }
        }

        size_t k = 0;
        float_t alpha = 1.0;
        base::DataVector dir(d);

        size_t breakIterationCounter = 0;
        const size_t BREAK_ITERATION_COUNTER_MAX = 10;

        while (k < N) {
          // calculate gradient
          fx = fGradient.eval(x, gradFx);
          k++;

          const float_t gradFxNorm = gradFx.l2Norm();

          if (gradFxNorm == 0.0) {
            break;
          }

          // DEBUG
          /*std::cout << "\nk = " << k << "\n";
          std::cout << "x = " << x.toString() << "\n";
          std::cout << "fx = " << fx << "\n";
          std::cout << "gradFx = " << gradFx.toString() << "\n";*/

          for (size_t i = 0; i < d; i++) {
            dir[i] = 0.0;

            for (size_t j = 0; j < d; j++) {
              dir[i] -= inverseHessian.get(i, j) * gradFx[j];
            }
          }

          if (dir.dotProduct(gradFx) > 0.0) {
            for (size_t t = 0; t < d; t++) {
              dir[t] = -gradFx[t] / gradFxNorm;
            }
          }

          for (size_t t = 0; t < d; t++) {
            xNew[t] = x[t] + alpha * dir[t];
          }

          // evaluate at new point
          fxNew = f.eval(xNew);
          k++;

          // inner product of gradient and search direction
          const float_t gradFxTimesDir = gradFx.dotProduct(dir);

          // line search
          while (fxNew > fx + rhoLs * alpha * gradFxTimesDir) {
            alpha *= rhoAlphaMinus;

            // recalculate new point
            for (size_t t = 0; t < d; t++) {
              xNew[t] = x[t] + alpha * dir[t];
            }

            // evaluate at new point
            fxNew = fGradient.eval(xNew, gradFxNew);
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

          // increase step size
          alpha *= rhoAlphaPlus;

          const float_t deltaTimesY = delta.dotProduct(y);
          // DEBUG
          /*std::cout << "xNew = " << xNew.toString() << "\n";
          std::cout << "fxNew = " << fxNew << "\n";
          std::cout << "gradFxNew = " << gradFxNew.toString() << "\n";
          std::cout << "delta = " << delta.toString() << "\n";
          std::cout << "y = " << y.toString() << "\n";
          std::cout << "deltaTimesY = " << deltaTimesY << "\n";*/

          if (deltaTimesY != 0.0) {
            for (size_t i = 0; i < d; i++) {
              for (size_t j = 0; j < d; j++) {
                M.set(i, j, (i == j ? 1.0 : 0.0) -
                      y[i] * delta[j] / deltaTimesY);
              }
            }

            for (size_t i = 0; i < d; i++) {
              for (size_t j = 0; j < d; j++) {
                float_t entry = delta[i] * delta[j] / deltaTimesY;

                for (size_t p = 0; p < d; p++) {
                  for (size_t q = 0; q < d; q++) {
                    entry +=
                      M.get(p, i) * inverseHessian.get(p, q) * M.get(q, j);
                  }
                }

                inverseHessianNew.set(i, j, entry);
              }
            }

            inverseHessian = inverseHessianNew;
          }

          // status printing
          printer.printStatusUpdate(
            std::to_string(k) + " evaluations, f(x) = " +
            std::to_string(fx));

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

        printer.printStatusUpdate(
          std::to_string(k) + " evaluations, f(x) = " +
          std::to_string(fx));;
        printer.printStatusEnd();

        return fx;
      }

      ObjectiveGradient& BFGS::getObjectiveGradient() const {
        return fGradient;
      }

      float_t BFGS::getTolerance() const {
        return theta;
      }

      void BFGS::setTolerance(float_t tolerance) {
        theta = tolerance;
      }

      float_t BFGS::getStepSizeIncreaseFactor() const {
        return rhoAlphaPlus;
      }

      void BFGS::setStepSizeIncreaseFactor(
        float_t stepSizeIncreaseFactor) {
        rhoAlphaPlus = stepSizeIncreaseFactor;
      }

      float_t BFGS::getStepSizeDecreaseFactor() const {
        return rhoAlphaMinus;
      }

      void BFGS::setStepSizeDecreaseFactor(
        float_t stepSizeDecreaseFactor) {
        rhoAlphaMinus = stepSizeDecreaseFactor;
      }

      float_t BFGS::getLineSearchAccuracy() const {
        return rhoLs;
      }

      void BFGS::setLineSearchAccuracy(
        float_t lineSearchAccuracy) {
        rhoLs = lineSearchAccuracy;
      }

    }
  }
}

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/constrained/LogBarrier.hpp>
#include <sgpp/optimization/function/EmptyConstraintFunction.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      namespace {
        class PenalizedObjectiveFunction : public ObjectiveFunction {
          public:
            PenalizedObjectiveFunction(ObjectiveFunction& f,
                                       ConstraintFunction& g,
                                       float_t mu) :
              ObjectiveFunction(f.getDimension()),
              f(f),
              g(g),
              mu(mu),
              m(g.getNumberOfConstraints()) {
            }

            float_t eval(const base::DataVector& x) {
              for (size_t t = 0; t < d; t++) {
                if ((x.get(t) < 0.0) || (x.get(t) > 1.0)) {
                  return INFINITY;
                }
              }

              const float_t fx = f.eval(x);
              base::DataVector gx(m);
              g.eval(x, gx);

              float_t value = fx;

              for (size_t i = 0; i < m; i++) {
                if (gx[i] < 0.0) {
                  value -= mu * std::log(-gx[i]);
                } else {
                  return INFINITY;
                }
              }

              return value;
            }

            void clone(std::unique_ptr<ObjectiveFunction>& clone) const {
              clone = std::unique_ptr<ObjectiveFunction>(
                        new PenalizedObjectiveFunction(*this));
            }

            float_t getMu() const {
              return mu;
            }

            void setMu(float_t mu) {
              this->mu = mu;
            }

          protected:
            ObjectiveFunction& f;
            ConstraintFunction& g;
            float_t mu;
            size_t m;
        };

        class PenalizedObjectiveGradient : public ObjectiveGradient {
          public:
            PenalizedObjectiveGradient(ObjectiveGradient& fGradient,
                                       ConstraintGradient& gGradient,
                                       float_t mu) :
              ObjectiveGradient(fGradient.getDimension()),
              fGradient(fGradient),
              gGradient(gGradient),
              mu(mu),
              m(gGradient.getNumberOfConstraints()) {
            }

            float_t eval(const base::DataVector& x,
                         base::DataVector& gradient) {
              for (size_t t = 0; t < d; t++) {
                if ((x.get(t) < 0.0) || (x.get(t) > 1.0)) {
                  return INFINITY;
                }
              }

              base::DataVector gradFx(d);
              const float_t fx = fGradient.eval(x, gradFx);
              base::DataVector gx(m);
              base::DataMatrix gradGx(m, d);
              gGradient.eval(x, gx, gradGx);

              float_t value = fx;
              gradient.resize(d);
              gradient = gradFx;

              for (size_t i = 0; i < m; i++) {
                const float_t gxi = gx[i];

                if (gxi < 0.0) {
                  value -= mu * std::log(-gx[i]);

                  for (size_t t = 0; t < d; t++) {
                    gradient[t] -= mu * gradGx.get(i, t) / gx[i];
                  }
                } else {
                  return INFINITY;
                }
              }

              return value;
            }

            void clone(std::unique_ptr<ObjectiveGradient>& clone) const {
              clone = std::unique_ptr<ObjectiveGradient>(
                        new PenalizedObjectiveGradient(*this));
            }

            float_t getMu() const {
              return mu;
            }

            void setMu(float_t mu) {
              this->mu = mu;
            }

          protected:
            ObjectiveGradient& fGradient;
            ConstraintGradient& gGradient;
            float_t mu;
            size_t m;
        };
      }

      LogBarrier::LogBarrier(
        ObjectiveFunction& f,
        ObjectiveGradient& fGradient,
        ConstraintFunction& g,
        ConstraintGradient& gGradient,
        size_t maxItCount,
        float_t tolerance,
        float_t barrierStartValue,
        float_t barrierDecreaseFactor) :
        ConstrainedOptimizer(f, g, emptyConstraintFunction, maxItCount),
        fGradient(fGradient),
        gGradient(gGradient),
        theta(tolerance),
        mu0(barrierStartValue),
        rhoMuMinus(barrierDecreaseFactor) {
      }

      float_t LogBarrier::optimize(base::DataVector& xOpt) {
        printer.printStatusBegin("Optimizing (Log Barrier)...");

        const size_t d = f.getDimension();

        base::DataVector x(x0);
        float_t fx = f.eval(x);

        base::DataVector xNew(d);

        float_t mu = mu0;

        size_t breakIterationCounter = 0;
        const size_t BREAK_ITERATION_COUNTER_MAX = 10;
        size_t k = 1;

        const size_t unconstrainedN = N / 20;

        PenalizedObjectiveFunction fPenalized(f, g, mu);
        PenalizedObjectiveGradient fPenalizedGradient(
          fGradient, gGradient, mu);

        // http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page3389.htm
        while (k < N) {
          fPenalized.setMu(mu);
          fPenalizedGradient.setMu(mu);

          AdaptiveGradientDescent unconstrainedOptimizer(
            fPenalized, fPenalizedGradient, unconstrainedN, 10.0 * theta);
          unconstrainedOptimizer.setStartingPoint(x);
          unconstrainedOptimizer.optimize(xNew);
          k += unconstrainedN;

          x = xNew;
          fx = f.eval(x);
          k++;

          // status printing
          printer.printStatusUpdate(
            std::to_string(k) + " evaluations, f(x) = " +
            std::to_string(fx));

          mu *= rhoMuMinus;

          xNew.sub(x);

          if (xNew.l2Norm() < theta) {
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
          std::to_string(fx));
        printer.printStatusEnd();

        return fx;
      }

      ObjectiveGradient& LogBarrier::getObjectiveGradient() const {
        return fGradient;
      }

      ConstraintGradient& LogBarrier::getInequalityConstraintGradient() const {
        return gGradient;
      }

      float_t LogBarrier::getTolerance() const {
        return theta;
      }

      void LogBarrier::setTolerance(float_t tolerance) {
        theta = tolerance;
      }

      float_t LogBarrier::getBarrierStartValue() const {
        return mu0;
      }

      void LogBarrier::setBarrierStartValue(float_t barrierStartValue) {
        mu0 = barrierStartValue;
      }

      float_t LogBarrier::getBarrierDecreaseFactor() const {
        return rhoMuMinus;
      }

      void LogBarrier::setBarrierDecreaseFactor(
        float_t barrierDecreaseFactor) {
        rhoMuMinus = barrierDecreaseFactor;
      }

    }
  }
}

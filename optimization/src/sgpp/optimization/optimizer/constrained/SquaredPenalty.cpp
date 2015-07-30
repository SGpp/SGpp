// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/constrained/SquaredPenalty.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      namespace {
        class PenalizedObjectiveFunction : public ObjectiveFunction {
          public:
            PenalizedObjectiveFunction(ObjectiveFunction& f,
                                       ConstraintFunction& g,
                                       ConstraintFunction& h,
                                       float_t mu) :
              ObjectiveFunction(f.getDimension()),
              f(f),
              g(g),
              h(h),
              mu(mu),
              mG(g.getNumberOfConstraints()),
              mH(h.getNumberOfConstraints()) {
            }

            float_t eval(const base::DataVector& x) {
              for (size_t t = 0; t < d; t++) {
                if ((x.get(t) < 0.0) || (x.get(t) > 1.0)) {
                  return INFINITY;
                }
              }

              const float_t fx = f.eval(x);

              base::DataVector gx(mG);
              g.eval(x, gx);

              base::DataVector hx(mH);
              h.eval(x, hx);

              // DEBUG
              /*std::cout << "x = " << x.toString() << "\n";
              std::cout << "fx = " << fx << "\n";
              std::cout << "gx = " << gx.toString() << "\n";
              std::cout << "hx = " << hx.toString() << "\n";*/

              float_t value = fx;

              for (size_t i = 0; i < mG; i++) {
                if (gx[i] > 0.0) {
                  value += mu * gx[i] * gx[i];
                }
              }

              for (size_t i = 0; i < mH; i++) {
                value += mu * hx[i] * hx[i];
              }

              // DEBUG
              //std::cout << "value = " << value << "\n";

              return value;
            }

            void clone(std::unique_ptr<ObjectiveFunction>& clone) const {
              clone = std::unique_ptr<ObjectiveFunction>(
                        new PenalizedObjectiveFunction(*this));
            }

            void setMu(float_t mu) {
              this->mu = mu;
            }

          protected:
            ObjectiveFunction& f;
            ConstraintFunction& g;
            ConstraintFunction& h;
            float_t mu;
            size_t mG;
            size_t mH;
        };

        class PenalizedObjectiveGradient : public ObjectiveGradient {
          public:
            PenalizedObjectiveGradient(ObjectiveGradient& fGradient,
                                       ConstraintGradient& gGradient,
                                       ConstraintGradient& hGradient,
                                       float_t mu) :
              ObjectiveGradient(fGradient.getDimension()),
              fGradient(fGradient),
              gGradient(gGradient),
              hGradient(hGradient),
              mu(mu),
              mG(gGradient.getNumberOfConstraints()),
              mH(hGradient.getNumberOfConstraints()) {
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

              base::DataVector gx(mG);
              base::DataMatrix gradGx(mG, d);
              gGradient.eval(x, gx, gradGx);

              base::DataVector hx(mH);
              base::DataMatrix gradHx(mH, d);
              hGradient.eval(x, hx, gradHx);

              // DEBUG
              /*std::cout << "mu = " << mu << "\n";
              std::cout << "mG = " << mG << "\n";
              std::cout << "mH = " << mH << "\n";
              std::cout << "x = " << x.toString() << "\n";
              std::cout << "fx = " << fx << "\n";
              std::cout << "gradFx = " << gradFx.toString() << "\n";
              std::cout << "gx = " << gx.toString() << "\n";
              std::cout << "gradGx = " << gradGx.toString() << "\n";
              std::cout << "hx = " << hx.toString() << "\n";
              std::cout << "gradHx = " << gradHx.toString() << "\n";*/

              float_t value = fx;
              gradient.resize(d);
              gradient = gradFx;

              for (size_t i = 0; i < mG; i++) {
                const float_t gxi = gx[i];

                if (gxi > 0.0) {
                  value += mu * gxi * gxi;

                  for (size_t t = 0; t < d; t++) {
                    gradient[t] += mu * 2.0 * gxi * gradGx.get(i, t);
                  }
                }
              }

              for (size_t i = 0; i < mH; i++) {
                const float_t hxi = hx[i];
                value += mu * hxi * hxi;

                for (size_t t = 0; t < d; t++) {
                  gradient[t] += mu * 2.0 * hxi * gradHx.get(i, t);
                }
              }

              // DEBUG
              /*std::cout << "value = " << value << "\n";
              std::cout << "gradient = " << gradient.toString() << "\n";*/

              return value;
            }

            void clone(std::unique_ptr<ObjectiveGradient>& clone) const {
              clone = std::unique_ptr<ObjectiveGradient>(
                        new PenalizedObjectiveGradient(*this));
            }

            void setMu(float_t mu) {
              this->mu = mu;
            }

          protected:
            ObjectiveGradient& fGradient;
            ConstraintGradient& gGradient;
            ConstraintGradient& hGradient;
            float_t mu;
            size_t mG;
            size_t mH;
        };
      }

      SquaredPenalty::SquaredPenalty(
        ObjectiveFunction& f,
        ObjectiveGradient& fGradient,
        ConstraintFunction& g,
        ConstraintGradient& gGradient,
        ConstraintFunction& h,
        ConstraintGradient& hGradient,
        size_t maxItCount,
        float_t xTolerance,
        float_t constraintTolerance,
        float_t penaltyStartValue,
        float_t penaltyIncreaseFactor) :
        ConstrainedOptimizer(f, g, h, maxItCount),
        fGradient(fGradient),
        gGradient(gGradient),
        hGradient(hGradient),
        theta(xTolerance),
        epsilon(constraintTolerance),
        mu0(penaltyStartValue),
        rhoMuPlus(penaltyIncreaseFactor) {
      }

      float_t SquaredPenalty::optimize(base::DataVector& xOpt) {
        printer.printStatusBegin("Optimizing (Squared Penalty)...");

        const size_t d = f.getDimension();
        const size_t mG = g.getNumberOfConstraints();
        const size_t mH = h.getNumberOfConstraints();

        base::DataVector x(x0);
        float_t fx = f.eval(x);

        base::DataVector xNew(d);

        base::DataVector gx(mG);
        base::DataVector hx(mH);

        float_t mu = mu0;

        size_t breakIterationCounter = 0;
        const size_t BREAK_ITERATION_COUNTER_MAX = 10;
        size_t k = 1;

        const size_t unconstrainedN = N / 20;

        PenalizedObjectiveFunction fPenalized(f, g, h, mu);
        PenalizedObjectiveGradient fPenalizedGradient(
          fGradient, gGradient, hGradient, mu);

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

          mu *= rhoMuPlus;

          g.eval(x, gx);
          h.eval(x, hx);

          xNew.sub(x);

          if ((xNew.l2Norm() < theta) &&
              (gx.max() < epsilon) &&
              (hx.maxNorm() < epsilon)) {
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

      ObjectiveGradient& SquaredPenalty::getObjectiveGradient() const {
        return fGradient;
      }

      ConstraintGradient&
      SquaredPenalty::getInequalityConstraintGradient() const {
        return gGradient;
      }

      ConstraintGradient& SquaredPenalty::getEqualityConstraintGradient() const {
        return hGradient;
      }

      float_t SquaredPenalty::getXTolerance() const {
        return theta;
      }

      void SquaredPenalty::setXTolerance(float_t xTolerance) {
        theta = xTolerance;
      }

      float_t SquaredPenalty::getConstraintTolerance() const {
        return epsilon;
      }

      void SquaredPenalty::setConstraintTolerance(float_t constraintTolerance) {
        epsilon = constraintTolerance;
      }

      float_t SquaredPenalty::getPenaltyStartValue() const {
        return mu0;
      }

      void SquaredPenalty::setPenaltyStartValue(float_t penaltyStartValue) {
        mu0 = penaltyStartValue;
      }

      float_t SquaredPenalty::getPenaltyIncreaseFactor() const {
        return rhoMuPlus;
      }

      void SquaredPenalty::setPenaltyIncreaseFactor(
        float_t penaltyIncreaseFactor) {
        rhoMuPlus = penaltyIncreaseFactor;
      }

    }
  }
}

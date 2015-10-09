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
        class PenalizedObjectiveFunction : public ScalarFunction {
          public:
            PenalizedObjectiveFunction(ScalarFunction& f,
                                       VectorFunction& g,
                                       VectorFunction& h,
                                       float_t mu) :
              ScalarFunction(f.getNumberOfParameters()),
              f(f),
              g(g),
              h(h),
              mu(mu),
              mG(g.getNumberOfComponents()),
              mH(h.getNumberOfComponents()) {
            }

            float_t eval(const base::DataVector& x) {
              for (size_t t = 0; t < d; t++) {
                if ((x[t] < 0.0) || (x[t] > 1.0)) {
                  return INFINITY;
                }
              }

              const float_t fx = f.eval(x);

              base::DataVector gx(mG);
              g.eval(x, gx);

              base::DataVector hx(mH);
              h.eval(x, hx);

              float_t value = fx;

              for (size_t i = 0; i < mG; i++) {
                if (gx[i] > 0.0) {
                  value += mu * gx[i] * gx[i];
                }
              }

              for (size_t i = 0; i < mH; i++) {
                value += mu * hx[i] * hx[i];
              }

              return value;
            }

            void clone(std::unique_ptr<ScalarFunction>& clone) const {
              clone = std::unique_ptr<ScalarFunction>(
                        new PenalizedObjectiveFunction(*this));
            }

            void setMu(float_t mu) {
              this->mu = mu;
            }

          protected:
            ScalarFunction& f;
            VectorFunction& g;
            VectorFunction& h;
            float_t mu;
            size_t mG;
            size_t mH;
        };

        class PenalizedObjectiveGradient : public ScalarFunctionGradient {
          public:
            PenalizedObjectiveGradient(ScalarFunctionGradient& fGradient,
                                       VectorFunctionGradient& gGradient,
                                       VectorFunctionGradient& hGradient,
                                       float_t mu) :
              ScalarFunctionGradient(fGradient.getNumberOfParameters()),
              fGradient(fGradient),
              gGradient(gGradient),
              hGradient(hGradient),
              mu(mu),
              mG(gGradient.getNumberOfComponents()),
              mH(hGradient.getNumberOfComponents()) {
            }

            float_t eval(const base::DataVector& x,
                         base::DataVector& gradient) {
              for (size_t t = 0; t < d; t++) {
                if ((x[t] < 0.0) || (x[t] > 1.0)) {
                  gradient.setAll(NAN);
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

              float_t value = fx;
              gradient.resize(d);
              gradient = gradFx;

              for (size_t i = 0; i < mG; i++) {
                const float_t gxi = gx[i];

                if (gxi > 0.0) {
                  value += mu * gxi * gxi;

                  for (size_t t = 0; t < d; t++) {
                    gradient[t] += mu * 2.0 * gxi * gradGx(i, t);
                  }
                }
              }

              for (size_t i = 0; i < mH; i++) {
                const float_t hxi = hx[i];
                value += mu * hxi * hxi;

                for (size_t t = 0; t < d; t++) {
                  gradient[t] += mu * 2.0 * hxi * gradHx(i, t);
                }
              }

              return value;
            }

            void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const {
              clone = std::unique_ptr<ScalarFunctionGradient>(
                        new PenalizedObjectiveGradient(*this));
            }

            void setMu(float_t mu) {
              this->mu = mu;
            }

          protected:
            ScalarFunctionGradient& fGradient;
            VectorFunctionGradient& gGradient;
            VectorFunctionGradient& hGradient;
            float_t mu;
            size_t mG;
            size_t mH;
        };
      }

      SquaredPenalty::SquaredPenalty(
        ScalarFunction& f,
        ScalarFunctionGradient& fGradient,
        VectorFunction& g,
        VectorFunctionGradient& gGradient,
        VectorFunction& h,
        VectorFunctionGradient& hGradient,
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
        rhoMuPlus(penaltyIncreaseFactor),
        kHist() {
      }

      void SquaredPenalty::optimize() {
        printer.printStatusBegin("Optimizing (Squared Penalty)...");

        const size_t d = f.getNumberOfParameters();

        xOpt.resize(0);
        fOpt = NAN;
        xHist.resize(0, d);
        fHist.resize(0);
        kHist.clear();

        const size_t mG = g.getNumberOfComponents();
        const size_t mH = h.getNumberOfComponents();

        base::DataVector x(x0);
        float_t fx = f.eval(x);

        xHist.appendRow(x);
        fHist.append(fx);

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
          unconstrainedOptimizer.optimize();
          xNew = unconstrainedOptimizer.getOptimalPoint();

          const size_t numberInnerEvaluations =
            unconstrainedOptimizer.getHistoryOfOptimalPoints().getNrows();
          k += numberInnerEvaluations;

          x = xNew;
          fx = f.eval(x);
          g.eval(x, gx);
          h.eval(x, hx);
          k++;

          xHist.appendRow(x);
          fHist.append(fx);
          kHist.push_back(numberInnerEvaluations);

          // status printing
          printer.printStatusUpdate(
            std::to_string(k) + " evaluations, x = " + x.toString() +
            ", f(x) = " + std::to_string(fx) +
            ", g(x) = " + gx.toString() +
            ", h(x) = " + hx.toString());

          mu *= rhoMuPlus;

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
        fOpt = fx;
        printer.printStatusEnd();
      }

      ScalarFunctionGradient& SquaredPenalty::getObjectiveGradient() const {
        return fGradient;
      }

      VectorFunctionGradient&
      SquaredPenalty::getInequalityConstraintGradient() const {
        return gGradient;
      }

      VectorFunctionGradient& SquaredPenalty::getEqualityConstraintGradient() const {
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

      const std::vector<size_t>&
      SquaredPenalty::getHistoryOfInnerIterations() const {
        return kHist;
      }

      void SquaredPenalty::clone(
        std::unique_ptr<UnconstrainedOptimizer>& clone) const {
        clone = std::unique_ptr<UnconstrainedOptimizer>(
                  new SquaredPenalty(*this));
      }
    }
  }
}

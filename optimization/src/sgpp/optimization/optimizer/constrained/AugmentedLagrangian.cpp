// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/constrained/AugmentedLagrangian.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>
#include <sgpp/optimization/function/vector/EmptyVectorFunction.hpp>
#include <sgpp/optimization/function/vector/EmptyVectorFunctionGradient.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      namespace {
        class PenalizedObjectiveFunction : public ScalarFunction {
          public:
            PenalizedObjectiveFunction(ScalarFunction& f,
                                       VectorFunction& g,
                                       VectorFunction& h,
                                       float_t mu,
                                       base::DataVector& lambda) :
              ScalarFunction(f.getNumberOfParameters()),
              f(f),
              g(g),
              h(h),
              mu(mu),
              lambda(lambda),
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
                if ((gx[i] >= 0.0) || (lambda[i] > 0)) {
                  value += (mu * gx[i] + lambda[i]) * gx[i];
                } else {
                  value += lambda[i] * gx[i];
                }
              }

              for (size_t i = 0; i < mH; i++) {
                value += (mu * hx[i] + lambda[mG + i]) * hx[i];
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
            base::DataVector& lambda;
            size_t mG;
            size_t mH;
        };

        class PenalizedObjectiveGradient : public ScalarFunctionGradient {
          public:
            PenalizedObjectiveGradient(ScalarFunctionGradient& fGradient,
                                       VectorFunctionGradient& gGradient,
                                       VectorFunctionGradient& hGradient,
                                       float_t mu,
                                       base::DataVector& lambda) :
              ScalarFunctionGradient(fGradient.getNumberOfParameters()),
              fGradient(fGradient),
              gGradient(gGradient),
              hGradient(hGradient),
              mu(mu),
              lambda(lambda),
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
                const float_t lambdai = lambda[i];

                if ((gx[i] >= 0.0) || (lambdai > 0)) {
                  value += (mu * gxi + lambdai) * gxi;

                  for (size_t t = 0; t < d; t++) {
                    gradient[t] += (mu * 2.0 * gxi + lambdai) * gradGx(i, t);
                  }
                } else {
                  value += lambdai * gxi;

                  for (size_t t = 0; t < d; t++) {
                    gradient[t] += lambdai * gradGx(i, t);
                  }
                }
              }

              for (size_t i = 0; i < mH; i++) {
                const float_t hxi = hx[i];
                const float_t lambdai = lambda[mG + i];

                value += (mu * hxi + lambdai) * hxi;

                for (size_t t = 0; t < d; t++) {
                  gradient[t] += (mu * 2.0 * hxi + lambdai) * gradHx(i, t);
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
            base::DataVector& lambda;
            size_t mG;
            size_t mH;
        };

        class AuxiliaryObjectiveFunction : public ScalarFunction {
          public:
            AuxiliaryObjectiveFunction(size_t d, float_t sMin, float_t sMax) :
              ScalarFunction(d + 1),
              sMin(sMin),
              sMax(sMax) {
            }

            float_t eval(const base::DataVector& x) {
              const size_t d = this->d - 1;

              for (size_t t = 0; t < d + 1; t++) {
                if ((x[t] < 0.0) || (x[t] > 1.0)) {
                  return INFINITY;
                }
              }

              return x[d] * (sMax - sMin) + sMin;
            }

            void clone(std::unique_ptr<ScalarFunction>& clone) const {
              clone = std::unique_ptr<ScalarFunction>(
                        new AuxiliaryObjectiveFunction(*this));
            }

          protected:
            float_t sMin;
            float_t sMax;
        };

        class AuxiliaryObjectiveGradient : public ScalarFunctionGradient {
          public:
            AuxiliaryObjectiveGradient(size_t d, float_t sMin, float_t sMax) :
              ScalarFunctionGradient(d + 1),
              sMin(sMin),
              sMax(sMax) {
            }

            float_t eval(const base::DataVector& x,
                         base::DataVector& gradient) {
              const size_t d = this->d - 1;

              for (size_t t = 0; t < d + 1; t++) {
                if ((x[t] < 0.0) || (x[t] > 1.0)) {
                  gradient.setAll(NAN);
                  return INFINITY;
                }

                gradient[t] = ((t < d) ? 0.0 : (sMax - sMin));
              }

              return x[d] * (sMax - sMin) + sMin;
            }

            void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const {
              clone = std::unique_ptr<ScalarFunctionGradient>(
                        new AuxiliaryObjectiveGradient(*this));
            }

          protected:
            float_t sMin;
            float_t sMax;
        };

        class AuxiliaryConstraintFunction : public VectorFunction {
          public:
            AuxiliaryConstraintFunction(size_t d,
                                        VectorFunction& g,
                                        VectorFunction& h,
                                        float_t sMin,
                                        float_t sMax) :
              VectorFunction(d + 1,
                             g.getNumberOfComponents() +
                             2 * h.getNumberOfComponents() + 1),
              g(g),
              h(h),
              mG(g.getNumberOfComponents()),
              mH(h.getNumberOfComponents()),
              sMin(sMin),
              sMax(sMax) {
            }

            void eval(const base::DataVector& x,
                      base::DataVector& value) {
              const size_t d = this->d - 1;
              base::DataVector xPart(d);

              for (size_t t = 0; t < d + 1; t++) {
                if ((x[t] < 0.0) || (x[t] > 1.0)) {
                  value.setAll(INFINITY);
                  return;
                }

                if (t < d) {
                  xPart[t] = x[t];
                }
              }

              const float_t s = x[d] * (sMax - sMin) + sMin;
              base::DataVector gx(mG), hx(mH);

              g.eval(xPart, gx);
              h.eval(xPart, hx);

              value[0] = -s;

              for (size_t i = 0; i < mG; i++) {
                value[i + 1] = gx[i] - s;
              }

              for (size_t i = 0; i < mH; i++) {
                value[2 * i + mG + 1] = hx[i] - s;
                value[2 * i + mG + 2] = -hx[i] - s;
              }
            }

            void clone(std::unique_ptr<VectorFunction>& clone) const {
              clone = std::unique_ptr<VectorFunction>(
                        new AuxiliaryConstraintFunction(*this));
            }

          protected:
            VectorFunction& g;
            VectorFunction& h;
            size_t mG;
            size_t mH;
            float_t sMin;
            float_t sMax;
        };

        class AuxiliaryConstraintGradient : public VectorFunctionGradient {
          public:
            AuxiliaryConstraintGradient(size_t d,
                                        VectorFunctionGradient& gGradient,
                                        VectorFunctionGradient& hGradient,
                                        float_t sMin,
                                        float_t sMax) :
              VectorFunctionGradient(d + 1,
                                     gGradient.getNumberOfComponents() +
                                     2 * hGradient.getNumberOfComponents() + 1),
              gGradient(gGradient),
              hGradient(hGradient),
              mG(gGradient.getNumberOfComponents()),
              mH(hGradient.getNumberOfComponents()),
              sMin(sMin),
              sMax(sMax) {
            }

            void eval(const base::DataVector& x,
                      base::DataVector& value,
                      base::DataMatrix& gradient) {
              const size_t d = this->d - 1;
              base::DataVector xPart(d);

              for (size_t t = 0; t < d + 1; t++) {
                if ((x[t] < 0.0) || (x[t] > 1.0)) {
                  value.setAll(INFINITY);
                  gradient.setAll(NAN);
                  return;
                }

                if (t < d) {
                  xPart[t] = x[t];
                }
              }

              const float_t s = x[d] * (sMax - sMin) + sMin;
              base::DataVector gx(mG), hx(mH);
              base::DataMatrix gxGradient(mG, d), hxGradient(mH, d);

              gGradient.eval(xPart, gx, gxGradient);
              hGradient.eval(xPart, hx, hxGradient);

              value[0] = -s;

              for (size_t t = 0; t < d; t++) {
                gradient(0, t) = 0.0;
              }

              gradient(0, d) = -(sMax - sMin);

              for (size_t i = 0; i < mG; i++) {
                value[i + 1] = gx[i] - s;

                for (size_t t = 0; t < d; t++) {
                  gradient(i + 1, t) = gxGradient(i, t);
                }

                gradient(i + 1, d) = -(sMax - sMin);
              }

              for (size_t i = 0; i < mH; i++) {
                const size_t j = 2 * i + mG + 1;
                value[j] = hx[i] - s;
                value[j + 1] = -hx[i] - s;

                for (size_t t = 0; t < d; t++) {
                  gradient(j, t) = hxGradient(i, t);
                  gradient(j + 1, t) = -hxGradient(i, t);
                }

                gradient(j, d) = -(sMax - sMin);
                gradient(j + 1, d) = -(sMax - sMin);
              }
            }

            void clone(std::unique_ptr<VectorFunctionGradient>& clone) const {
              clone = std::unique_ptr<VectorFunctionGradient>(
                        new AuxiliaryConstraintGradient(*this));
            }

          protected:
            VectorFunctionGradient& gGradient;
            VectorFunctionGradient& hGradient;
            size_t mG;
            size_t mH;
            float_t sMin;
            float_t sMax;
        };
      }

      AugmentedLagrangian::AugmentedLagrangian(
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

      void AugmentedLagrangian::optimize() {
        printer.printStatusBegin("Optimizing (Augmented Lagrangian)...");

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
        base::DataVector lambda(mG + mH, 0.0);

        size_t breakIterationCounter = 0;
        const size_t BREAK_ITERATION_COUNTER_MAX = 10;
        size_t k = 1;

        const size_t unconstrainedN = N / 20;

        PenalizedObjectiveFunction fPenalized(f, g, h, mu, lambda);
        PenalizedObjectiveGradient fPenalizedGradient(
          fGradient, gGradient, hGradient, mu, lambda);

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

          for (size_t i = 0; i < mG; i++) {
            lambda[i] = std::max(lambda[i] + 2.0 * mu * gx[i], 0.0);
          }

          for (size_t i = 0; i < mH; i++) {
            lambda[mG + i] += 2.0 * mu * hx[i];
          }

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

      base::DataVector AugmentedLagrangian::findFeasiblePoint() const {
        const size_t d = f.getNumberOfParameters();
        const size_t mG = g.getNumberOfComponents();
        const size_t mH = h.getNumberOfComponents();
        base::DataVector x(d, 0.5);
        base::DataVector gx(mG);
        base::DataVector hx(mH);

        g.eval(x, gx);
        h.eval(x, hx);
        const float_t s0 = 1.1 * std::max(gx.max(), hx.maxNorm());

        if (s0 == 0.0) {
          return x;
        }

        const float_t sMin = -0.1 * s0;
        const float_t sMax = 1.1 * s0;

        AuxiliaryObjectiveFunction auxObjFun(d, sMin, sMax);
        AuxiliaryObjectiveGradient auxObjGrad(d, sMin, sMax);
        AuxiliaryConstraintFunction auxConstrFun(d, g, h, sMin, sMax);
        AuxiliaryConstraintGradient auxConstrGrad(
          d, gGradient, hGradient, sMin, sMax);

        base::DataVector auxX(d + 1);

        for (size_t t = 0; t < d; t++) {
          auxX[t] = x[t];
        }

        auxX[d] = (s0 - sMin) / (sMax - sMin);

        AugmentedLagrangian optimizer(
          auxObjFun, auxObjGrad, auxConstrFun, auxConstrGrad,
          emptyVectorFunction, emptyVectorFunctionGradient);
        optimizer.setStartingPoint(auxX);
        optimizer.optimize();
        auxX = optimizer.getOptimalPoint();

        for (size_t t = 0; t < d; t++) {
          x[t] = auxX[t];
        }

        return x;
      }

      ScalarFunctionGradient& AugmentedLagrangian::getObjectiveGradient() const {
        return fGradient;
      }

      VectorFunctionGradient&
      AugmentedLagrangian::getInequalityConstraintGradient() const {
        return gGradient;
      }

      VectorFunctionGradient&
      AugmentedLagrangian::getEqualityConstraintGradient() const {
        return hGradient;
      }

      float_t AugmentedLagrangian::getXTolerance() const {
        return theta;
      }

      void AugmentedLagrangian::setXTolerance(float_t xTolerance) {
        theta = xTolerance;
      }

      float_t AugmentedLagrangian::getConstraintTolerance() const {
        return epsilon;
      }

      void AugmentedLagrangian::setConstraintTolerance(
        float_t constraintTolerance) {
        epsilon = constraintTolerance;
      }

      float_t AugmentedLagrangian::getPenaltyStartValue() const {
        return mu0;
      }

      void AugmentedLagrangian::setPenaltyStartValue(
        float_t penaltyStartValue) {
        mu0 = penaltyStartValue;
      }

      float_t AugmentedLagrangian::getPenaltyIncreaseFactor() const {
        return rhoMuPlus;
      }

      void AugmentedLagrangian::setPenaltyIncreaseFactor(
        float_t penaltyIncreaseFactor) {
        rhoMuPlus = penaltyIncreaseFactor;
      }

      const std::vector<size_t>&
      AugmentedLagrangian::getHistoryOfInnerIterations() const {
        return kHist;
      }

      void AugmentedLagrangian::clone(
        std::unique_ptr<UnconstrainedOptimizer>& clone) const {
        clone = std::unique_ptr<UnconstrainedOptimizer>(
                  new AugmentedLagrangian(*this));
      }
    }
  }
}

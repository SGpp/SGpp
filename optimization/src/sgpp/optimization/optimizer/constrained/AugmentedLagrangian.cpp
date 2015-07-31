// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/constrained/AugmentedLagrangian.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>
#include <sgpp/optimization/function/EmptyConstraintFunction.hpp>
#include <sgpp/optimization/function/EmptyConstraintGradient.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      namespace {
        class PenalizedObjectiveFunction : public ObjectiveFunction {
          public:
            PenalizedObjectiveFunction(ObjectiveFunction& f,
                                       ConstraintFunction& g,
                                       ConstraintFunction& h,
                                       float_t mu,
                                       base::DataVector& lambda) :
              ObjectiveFunction(f.getDimension()),
              f(f),
              g(g),
              h(h),
              mu(mu),
              lambda(lambda),
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

              float_t value = fx;

              for (size_t i = 0; i < mG; i++) {
                if ((gx[i] >= 0.0) || (lambda[i] > 0)) {
                  value += (mu * gx[i] + lambda[i]) * gx[i];
                } else {
                  value += lambda[i] * gx[i];
                }
              }

              for (size_t i = 0; i < mH; i++) {
                value += (mu * hx[i] + lambda[i]) * hx[i];
              }

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
            base::DataVector& lambda;
            size_t mG;
            size_t mH;
        };

        class PenalizedObjectiveGradient : public ObjectiveGradient {
          public:
            PenalizedObjectiveGradient(ObjectiveGradient& fGradient,
                                       ConstraintGradient& gGradient,
                                       ConstraintGradient& hGradient,
                                       float_t mu,
                                       base::DataVector& lambda) :
              ObjectiveGradient(fGradient.getDimension()),
              fGradient(fGradient),
              gGradient(gGradient),
              hGradient(hGradient),
              mu(mu),
              lambda(lambda),
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

              float_t value = fx;
              gradient.resize(d);
              gradient = gradFx;

              for (size_t i = 0; i < mG; i++) {
                const float_t gxi = gx[i];
                const float_t lambdai = lambda[i];

                if ((gx[i] >= 0.0) || (lambdai > 0)) {
                  value += (mu * gxi + lambdai) * gxi;

                  for (size_t t = 0; t < d; t++) {
                    gradient[t] += (mu * 2.0 * gxi + lambdai) *
                                   gradGx.get(i, t);
                  }
                } else {
                  value += lambdai * gxi;

                  for (size_t t = 0; t < d; t++) {
                    gradient[t] += lambdai * gradGx.get(i, t);
                  }
                }
              }

              for (size_t i = 0; i < mH; i++) {
                const float_t hxi = hx[i];
                const float_t lambdai = lambda[i];

                value += (mu * hxi + lambdai) * hxi;

                for (size_t t = 0; t < d; t++) {
                  gradient[t] += (mu * 2.0 * hxi + lambdai) *
                                 gradHx.get(i, t);
                }
              }

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
            base::DataVector& lambda;
            size_t mG;
            size_t mH;
        };

        class AuxiliaryObjectiveFunction : public ObjectiveFunction {
          public:
            AuxiliaryObjectiveFunction(size_t d, float_t sMin, float_t sMax) :
              ObjectiveFunction(d + 1),
              sMin(sMin),
              sMax(sMax) {
            }

            float_t eval(const base::DataVector& x) {
              const size_t d = this->d - 1;

              for (size_t t = 0; t < d + 1; t++) {
                if ((x.get(t) < 0.0) || (x.get(t) > 1.0)) {
                  return INFINITY;
                }
              }

              return x.get(d) * (sMax - sMin) + sMin;
            }

            void clone(std::unique_ptr<ObjectiveFunction>& clone) const {
              clone = std::unique_ptr<ObjectiveFunction>(
                        new AuxiliaryObjectiveFunction(*this));
            }

          protected:
            float_t sMin;
            float_t sMax;
        };

        class AuxiliaryObjectiveGradient : public ObjectiveGradient {
          public:
            AuxiliaryObjectiveGradient(size_t d, float_t sMin, float_t sMax) :
              ObjectiveGradient(d + 1),
              sMin(sMin),
              sMax(sMax) {
            }

            float_t eval(const base::DataVector& x,
                         base::DataVector& gradient) {
              const size_t d = this->d - 1;

              for (size_t t = 0; t < d + 1; t++) {
                if ((x.get(t) < 0.0) || (x.get(t) > 1.0)) {
                  return INFINITY;
                }

                gradient[t] = ((t < d) ? 0.0 : (sMax - sMin));
              }

              return x.get(d) * (sMax - sMin) + sMin;
            }

            void clone(std::unique_ptr<ObjectiveGradient>& clone) const {
              clone = std::unique_ptr<ObjectiveGradient>(
                        new AuxiliaryObjectiveGradient(*this));
            }

          protected:
            float_t sMin;
            float_t sMax;
        };

        class AuxiliaryConstraintFunction : public ConstraintFunction {
          public:
            AuxiliaryConstraintFunction(size_t d,
                                        ConstraintFunction& g,
                                        ConstraintFunction& h,
                                        float_t sMin,
                                        float_t sMax) :
              ConstraintFunction(d + 1,
                                 g.getNumberOfConstraints() +
                                 2 * h.getNumberOfConstraints() + 1),
              g(g),
              h(h),
              mG(g.getNumberOfConstraints()),
              mH(h.getNumberOfConstraints()),
              sMin(sMin),
              sMax(sMax) {
            }

            void eval(const base::DataVector& x,
                      base::DataVector& value) {
              const size_t d = this->d - 1;
              base::DataVector xPart(d);

              for (size_t t = 0; t < d + 1; t++) {
                if ((x.get(t) < 0.0) || (x.get(t) > 1.0)) {
                  return;
                }

                if (t < d) {
                  xPart[t] = x.get(t);
                }
              }

              const float_t s = x.get(d) * (sMax - sMin) + sMin;
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

            void clone(std::unique_ptr<ConstraintFunction>& clone) const {
              clone = std::unique_ptr<ConstraintFunction>(
                        new AuxiliaryConstraintFunction(*this));
            }

          protected:
            ConstraintFunction& g;
            ConstraintFunction& h;
            size_t mG;
            size_t mH;
            float_t sMin;
            float_t sMax;
        };

        class AuxiliaryConstraintGradient : public ConstraintGradient {
          public:
            AuxiliaryConstraintGradient(size_t d,
                                        ConstraintGradient& gGradient,
                                        ConstraintGradient& hGradient,
                                        float_t sMin,
                                        float_t sMax) :
              ConstraintGradient(d + 1,
                                 gGradient.getNumberOfConstraints() +
                                 2 * hGradient.getNumberOfConstraints() + 1),
              gGradient(gGradient),
              hGradient(hGradient),
              mG(gGradient.getNumberOfConstraints()),
              mH(hGradient.getNumberOfConstraints()),
              sMin(sMin),
              sMax(sMax) {
            }

            void eval(const base::DataVector& x,
                      base::DataVector& value,
                      base::DataMatrix& gradient) {
              const size_t d = this->d - 1;
              base::DataVector xPart(d);
              //std::cout << "AuxiliaryConstraintGradient.eval(): x = " << x.toString() << ", d = " << d << "\n";

              for (size_t t = 0; t < d + 1; t++) {
                if ((x.get(t) < 0.0) || (x.get(t) > 1.0)) {
                  return;
                }

                if (t < d) {
                  xPart[t] = x.get(t);
                }
              }

              const float_t s = x.get(d) * (sMax - sMin) + sMin;
              base::DataVector gx(mG), hx(mH);
              base::DataMatrix gxGradient(mG, d), hxGradient(mH, d);

              gGradient.eval(xPart, gx, gxGradient);
              hGradient.eval(xPart, hx, hxGradient);

              value[0] = -s;

              for (size_t t = 0; t < d; t++) {
                gradient.set(0, t, 0.0);
              }

              gradient.set(0, d, -(sMax - sMin));

              for (size_t i = 0; i < mG; i++) {
                value[i + 1] = gx[i] - s;

                for (size_t t = 0; t < d; t++) {
                  gradient.set(i + 1, t, gxGradient.get(i, t));
                }

                gradient.set(i + 1, d, -(sMax - sMin));
              }

              for (size_t i = 0; i < mH; i++) {
                const size_t j = 2 * i + mG + 1;
                value[j] = hx[i] - s;
                value[j + 1] = -hx[i] - s;

                for (size_t t = 0; t < d; t++) {
                  gradient.set(j, t, hxGradient.get(i, t));
                  gradient.set(j + 1, t, -hxGradient.get(i, t));
                }

                gradient.set(j, d, -(sMax - sMin));
                gradient.set(j + 1, d, -(sMax - sMin));
              }
            }

            void clone(std::unique_ptr<ConstraintGradient>& clone) const {
              clone = std::unique_ptr<ConstraintGradient>(
                        new AuxiliaryConstraintGradient(*this));
            }

          protected:
            ConstraintGradient& gGradient;
            ConstraintGradient& hGradient;
            size_t mG;
            size_t mH;
            float_t sMin;
            float_t sMax;
        };
      }

      AugmentedLagrangian::AugmentedLagrangian(
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

      float_t AugmentedLagrangian::optimize(base::DataVector& xOpt) {
        printer.printStatusBegin("Optimizing (Augmented Lagrangian)...");

        const size_t d = f.getDimension();
        const size_t mG = g.getNumberOfConstraints();
        const size_t mH = h.getNumberOfConstraints();

        base::DataVector x(x0);
        float_t fx = f.eval(x);

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
          unconstrainedOptimizer.optimize(xNew);
          k += unconstrainedN;

          x = xNew;
          fx = f.eval(x);
          k++;

          // status printing
          printer.printStatusUpdate(
            std::to_string(k) + " evaluations, f(x) = " +
            std::to_string(fx));

          g.eval(x, gx);
          h.eval(x, hx);

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

        printer.printStatusUpdate(
          std::to_string(k) + " evaluations, f(x) = " +
          std::to_string(fx));
        printer.printStatusEnd();

        return fx;
      }

      base::DataVector AugmentedLagrangian::findFeasiblePoint() const {
        const size_t d = f.getDimension();
        const size_t mG = g.getNumberOfConstraints();
        const size_t mH = h.getNumberOfConstraints();
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
          emptyConstraintFunction, emptyConstraintGradient);
        optimizer.setStartingPoint(auxX);
        optimizer.optimize(auxX);

        for (size_t t = 0; t < d; t++) {
          x[t] = auxX[t];
        }

        return x;
      }

      ObjectiveGradient& AugmentedLagrangian::getObjectiveGradient() const {
        return fGradient;
      }

      ConstraintGradient&
      AugmentedLagrangian::getInequalityConstraintGradient() const {
        return gGradient;
      }

      ConstraintGradient&
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

    }
  }
}

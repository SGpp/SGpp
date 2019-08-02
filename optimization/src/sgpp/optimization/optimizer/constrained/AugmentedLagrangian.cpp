// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/function/vector/EmptyVectorFunction.hpp>
#include <sgpp/base/function/vector/EmptyVectorFunctionGradient.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/constrained/AugmentedLagrangian.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>

#include <algorithm>
#include <vector>

namespace sgpp {
namespace optimization {
namespace optimizer {

namespace {
class PenalizedObjectiveFunction : public base::ScalarFunction {
 public:
  PenalizedObjectiveFunction(base::ScalarFunction& f, base::VectorFunction& g,
                             base::VectorFunction& h, double mu, base::DataVector& lambda)
      : base::ScalarFunction(f.getNumberOfParameters()),
        f(f),
        g(g),
        h(h),
        mu(mu),
        lambda(lambda),
        mG(g.getNumberOfComponents()),
        mH(h.getNumberOfComponents()) {}

  double eval(const base::DataVector& x) {
    for (size_t t = 0; t < d; t++) {
      if ((x[t] < 0.0) || (x[t] > 1.0)) {
        return INFINITY;
      }
    }

    const double fx = f.eval(x);

    base::DataVector gx(mG);
    g.eval(x, gx);

    base::DataVector hx(mH);
    h.eval(x, hx);

    double value = fx;

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

  void clone(std::unique_ptr<base::ScalarFunction>& clone) const {
    clone = std::unique_ptr<base::ScalarFunction>(new PenalizedObjectiveFunction(*this));
  }

  void setMu(double mu) { this->mu = mu; }

 protected:
  base::ScalarFunction& f;
  base::VectorFunction& g;
  base::VectorFunction& h;
  double mu;
  base::DataVector& lambda;
  size_t mG;
  size_t mH;
};

class PenalizedObjectiveGradient : public base::ScalarFunctionGradient {
 public:
  PenalizedObjectiveGradient(base::ScalarFunctionGradient& fGradient,
                             base::VectorFunctionGradient& gGradient,
                             base::VectorFunctionGradient& hGradient, double mu,
                             base::DataVector& lambda)
      : base::ScalarFunctionGradient(fGradient.getNumberOfParameters()),
        fGradient(fGradient),
        gGradient(gGradient),
        hGradient(hGradient),
        mu(mu),
        lambda(lambda),
        mG(gGradient.getNumberOfComponents()),
        mH(hGradient.getNumberOfComponents()) {}

  double eval(const base::DataVector& x, base::DataVector& gradient) {
    for (size_t t = 0; t < d; t++) {
      if ((x[t] < 0.0) || (x[t] > 1.0)) {
        gradient.setAll(NAN);
        return INFINITY;
      }
    }

    base::DataVector gradFx(d);
    const double fx = fGradient.eval(x, gradFx);

    base::DataVector gx(mG);
    base::DataMatrix gradGx(mG, d);
    gGradient.eval(x, gx, gradGx);

    base::DataVector hx(mH);
    base::DataMatrix gradHx(mH, d);
    hGradient.eval(x, hx, gradHx);

    double value = fx;
    gradient.resize(d);
    gradient = gradFx;

    for (size_t i = 0; i < mG; i++) {
      const double gxi = gx[i];
      const double lambdai = lambda[i];

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
      const double hxi = hx[i];
      const double lambdai = lambda[mG + i];

      value += (mu * hxi + lambdai) * hxi;

      for (size_t t = 0; t < d; t++) {
        gradient[t] += (mu * 2.0 * hxi + lambdai) * gradHx(i, t);
      }
    }

    return value;
  }

  void clone(std::unique_ptr<base::ScalarFunctionGradient>& clone) const {
    clone = std::unique_ptr<base::ScalarFunctionGradient>(new PenalizedObjectiveGradient(*this));
  }

  void setMu(double mu) { this->mu = mu; }

 protected:
  base::ScalarFunctionGradient& fGradient;
  base::VectorFunctionGradient& gGradient;
  base::VectorFunctionGradient& hGradient;
  double mu;
  base::DataVector& lambda;
  size_t mG;
  size_t mH;
};

class AuxiliaryObjectiveFunction : public base::ScalarFunction {
 public:
  AuxiliaryObjectiveFunction(size_t d, double sMin, double sMax)
      : base::ScalarFunction(d + 1), sMin(sMin), sMax(sMax) {}

  double eval(const base::DataVector& x) {
    const size_t d = this->d - 1;

    for (size_t t = 0; t < d + 1; t++) {
      if ((x[t] < 0.0) || (x[t] > 1.0)) {
        return INFINITY;
      }
    }

    return x[d] * (sMax - sMin) + sMin;
  }

  void clone(std::unique_ptr<base::ScalarFunction>& clone) const {
    clone = std::unique_ptr<base::ScalarFunction>(new AuxiliaryObjectiveFunction(*this));
  }

 protected:
  double sMin;
  double sMax;
};

class AuxiliaryObjectiveGradient : public base::ScalarFunctionGradient {
 public:
  AuxiliaryObjectiveGradient(size_t d, double sMin, double sMax)
      : base::ScalarFunctionGradient(d + 1), sMin(sMin), sMax(sMax) {}

  double eval(const base::DataVector& x, base::DataVector& gradient) {
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

  void clone(std::unique_ptr<base::ScalarFunctionGradient>& clone) const {
    clone = std::unique_ptr<base::ScalarFunctionGradient>(new AuxiliaryObjectiveGradient(*this));
  }

 protected:
  double sMin;
  double sMax;
};

class AuxiliaryConstraintFunction : public base::VectorFunction {
 public:
  AuxiliaryConstraintFunction(size_t d, base::VectorFunction& g, base::VectorFunction& h,
                              double sMin, double sMax)
      : base::VectorFunction(d + 1, g.getNumberOfComponents() + 2 * h.getNumberOfComponents() + 1),
        g(g),
        h(h),
        mG(g.getNumberOfComponents()),
        mH(h.getNumberOfComponents()),
        sMin(sMin),
        sMax(sMax) {}

  void eval(const base::DataVector& x, base::DataVector& value) {
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

    const double s = x[d] * (sMax - sMin) + sMin;
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

  void clone(std::unique_ptr<base::VectorFunction>& clone) const {
    clone = std::unique_ptr<base::VectorFunction>(new AuxiliaryConstraintFunction(*this));
  }

 protected:
  base::VectorFunction& g;
  base::VectorFunction& h;
  size_t mG;
  size_t mH;
  double sMin;
  double sMax;
};

class AuxiliaryConstraintGradient : public base::VectorFunctionGradient {
 public:
  AuxiliaryConstraintGradient(size_t d, base::VectorFunctionGradient& gGradient,
                              base::VectorFunctionGradient& hGradient, double sMin, double sMax)
      : base::VectorFunctionGradient(
            d + 1, gGradient.getNumberOfComponents() + 2 * hGradient.getNumberOfComponents() + 1),
        gGradient(gGradient),
        hGradient(hGradient),
        mG(gGradient.getNumberOfComponents()),
        mH(hGradient.getNumberOfComponents()),
        sMin(sMin),
        sMax(sMax) {}

  void eval(const base::DataVector& x, base::DataVector& value, base::DataMatrix& gradient) {
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

    const double s = x[d] * (sMax - sMin) + sMin;
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

  void clone(std::unique_ptr<base::VectorFunctionGradient>& clone) const {
    clone = std::unique_ptr<base::VectorFunctionGradient>(new AuxiliaryConstraintGradient(*this));
  }

 protected:
  base::VectorFunctionGradient& gGradient;
  base::VectorFunctionGradient& hGradient;
  size_t mG;
  size_t mH;
  double sMin;
  double sMax;
};
}  // namespace

AugmentedLagrangian::AugmentedLagrangian(
    const base::ScalarFunction& f,
    const base::VectorFunction& g,
    const base::VectorFunction& h, size_t maxItCount,
    double xTolerance, double constraintTolerance, double penaltyStartValue,
    double penaltyIncreaseFactor)
    : ConstrainedOptimizer(f, g, h, maxItCount),
      theta(xTolerance),
      epsilon(constraintTolerance),
      mu0(penaltyStartValue),
      rhoMuPlus(penaltyIncreaseFactor),
      xHistInner(0, 0),
      kHistInner() {
}

AugmentedLagrangian::AugmentedLagrangian(
    const base::ScalarFunction& f, const base::ScalarFunctionGradient& fGradient,
    const base::VectorFunction& g, const base::VectorFunctionGradient& gGradient,
    const base::VectorFunction& h, const base::VectorFunctionGradient& hGradient, size_t maxItCount,
    double xTolerance, double constraintTolerance, double penaltyStartValue,
    double penaltyIncreaseFactor)
    : ConstrainedOptimizer(f, fGradient, g, gGradient, h, hGradient, maxItCount),
      theta(xTolerance),
      epsilon(constraintTolerance),
      mu0(penaltyStartValue),
      rhoMuPlus(penaltyIncreaseFactor),
      xHistInner(0, 0),
      kHistInner() {
  dynamic_cast<AdaptiveGradientDescent*>(
      unconstrainedOptimizer.get())->setTolerance(10.0 * theta);
}

AugmentedLagrangian::AugmentedLagrangian(
    const UnconstrainedOptimizer& unconstrainedOptimizer,
    const base::VectorFunction& g,
    const base::VectorFunctionGradient* gGradient,
    const base::VectorFunction& h,
    const base::VectorFunctionGradient* hGradient, size_t maxItCount,
    double xTolerance, double constraintTolerance, double penaltyStartValue,
    double penaltyIncreaseFactor)
    : ConstrainedOptimizer(unconstrainedOptimizer, g, gGradient, h, hGradient, maxItCount),
      theta(xTolerance),
      epsilon(constraintTolerance),
      mu0(penaltyStartValue),
      rhoMuPlus(penaltyIncreaseFactor),
      xHistInner(0, 0),
      kHistInner() {
}

AugmentedLagrangian::AugmentedLagrangian(const AugmentedLagrangian& other)
    : ConstrainedOptimizer(other),
      theta(other.theta),
      epsilon(other.epsilon),
      mu0(other.mu0),
      rhoMuPlus(other.rhoMuPlus),
      xHistInner(other.xHistInner),
      kHistInner(other.kHistInner) {
}

AugmentedLagrangian::~AugmentedLagrangian() {}

void AugmentedLagrangian::optimize() {
  base::Printer::getInstance().printStatusBegin("Optimizing (Augmented Lagrangian)...");

  const size_t d = f->getNumberOfParameters();

  xOpt.resize(0);
  fOpt = NAN;
  xHist.resize(0, d);
  fHist.resize(0);
  xHistInner.resize(0, d);
  kHistInner.clear();

  const size_t mG = g->getNumberOfComponents();
  const size_t mH = h->getNumberOfComponents();

  base::DataVector x(x0);
  double fx = f->eval(x);

  xHist.appendRow(x);
  fHist.append(fx);

  base::DataVector xOld(d);

  base::DataVector gx(mG);
  base::DataVector hx(mH);

  double mu = mu0;
  base::DataVector lambda(mG + mH, 0.0);

  size_t breakIterationCounter = 0;
  const size_t BREAK_ITERATION_COUNTER_MAX = 10;
  size_t k = 1;

  const size_t unconstrainedN = N / 20;

  PenalizedObjectiveFunction fPenalized(*f, *g, *h, mu, lambda);
  std::unique_ptr<PenalizedObjectiveGradient> fPenalizedGradient;

  if ((fGradient != nullptr) && (gGradient != nullptr) && (hGradient != nullptr)) {
    fPenalizedGradient.reset(
        new PenalizedObjectiveGradient(*fGradient, *gGradient, *hGradient, mu, lambda));
  }

  while (k < N) {
    fPenalized.setMu(mu);
    unconstrainedOptimizer->setObjectiveFunction(fPenalized);

    if (fPenalizedGradient) {
      fPenalizedGradient->setMu(mu);
      unconstrainedOptimizer->setObjectiveGradient(fPenalizedGradient.get());
    }

    unconstrainedOptimizer->setN(unconstrainedN);
    unconstrainedOptimizer->setStartingPoint(x);
    unconstrainedOptimizer->optimize();
    x = unconstrainedOptimizer->getOptimalPoint();

    const base::DataMatrix& innerPoints = unconstrainedOptimizer->getHistoryOfOptimalPoints();
    const size_t numberInnerIterations = innerPoints.getNrows();
    k += numberInnerIterations;

    fx = f->eval(x);
    g->eval(x, gx);
    h->eval(x, hx);
    k++;

    xHist.appendRow(x);
    fHist.append(fx);
    xHistInner.resize(xHistInner.getNrows() + numberInnerIterations, d);

    for (size_t i = 0; i < numberInnerIterations; i++) {
      for (size_t t = 0; t < d; t++) {
        xHistInner(xHistInner.getNrows() - numberInnerIterations + i, t) = innerPoints(i, t);
      }
    }

    kHistInner.push_back(numberInnerIterations);

    // status printing
    base::Printer::getInstance().printStatusUpdate(
        std::to_string(k) + " evaluations, x = " + x.toString() + ", f(x) = " + std::to_string(fx) +
        ", g(x) = " + gx.toString() + ", h(x) = " + hx.toString());

    for (size_t i = 0; i < mG; i++) {
      lambda[i] = std::max(lambda[i] + 2.0 * mu * gx[i], 0.0);
    }

    for (size_t i = 0; i < mH; i++) {
      lambda[mG + i] += 2.0 * mu * hx[i];
    }

    mu *= rhoMuPlus;

    xOld.sub(x);

    if ((xOld.l2Norm() < theta) && (gx.max() < epsilon) && (hx.maxNorm() < epsilon)) {
      breakIterationCounter++;

      if (breakIterationCounter >= BREAK_ITERATION_COUNTER_MAX) {
        break;
      }
    } else {
      breakIterationCounter = 0;
    }

    xOld = x;
  }

  xOpt.resize(d);
  xOpt = x;
  fOpt = fx;
  base::Printer::getInstance().printStatusEnd();
}

base::DataVector AugmentedLagrangian::findFeasiblePoint() const {
  const size_t d = f->getNumberOfParameters();
  const size_t mG = g->getNumberOfComponents();
  const size_t mH = h->getNumberOfComponents();
  base::DataVector x(d, 0.5);
  base::DataVector gx(mG);
  base::DataVector hx(mH);

  g->eval(x, gx);
  h->eval(x, hx);
  const double s0 = 1.1 * std::max(gx.max(), hx.maxNorm());

  if (s0 == 0.0) {
    return x;
  }

  const double sMin = -0.1 * s0;
  const double sMax = 1.1 * s0;

  AuxiliaryObjectiveFunction auxObjFun(d, sMin, sMax);
  AuxiliaryConstraintFunction auxConstrFun(d, *g, *h, sMin, sMax);
  std::unique_ptr<AuxiliaryObjectiveGradient> auxObjGrad;
  std::unique_ptr<AuxiliaryConstraintGradient> auxConstrGrad;

  if ((gGradient != nullptr) && (hGradient != nullptr)) {
    auxObjGrad.reset(new AuxiliaryObjectiveGradient(d, sMin, sMax));
    auxConstrGrad.reset(new AuxiliaryConstraintGradient(
        d, *gGradient, *hGradient, sMin, sMax));
  }

  base::DataVector auxX(d + 1);

  for (size_t t = 0; t < d; t++) {
    auxX[t] = x[t];
  }

  auxX[d] = (s0 - sMin) / (sMax - sMin);

  std::unique_ptr<AugmentedLagrangian> optimizer;

  if (auxObjGrad) {
    optimizer.reset(new AugmentedLagrangian(
        auxObjFun, *auxObjGrad, auxConstrFun, *auxConstrGrad,
        base::EmptyVectorFunction::getInstance(),
        base::EmptyVectorFunctionGradient::getInstance()));
  } else {
    optimizer.reset(new AugmentedLagrangian(
        auxObjFun, auxConstrFun, base::EmptyVectorFunction::getInstance()));
  }

  optimizer->setStartingPoint(auxX);
  optimizer->optimize();
  auxX = optimizer->getOptimalPoint();

  for (size_t t = 0; t < d; t++) {
    x[t] = auxX[t];
  }

  return x;
}

double AugmentedLagrangian::getXTolerance() const { return theta; }

void AugmentedLagrangian::setXTolerance(double xTolerance) { theta = xTolerance; }

double AugmentedLagrangian::getConstraintTolerance() const { return epsilon; }

void AugmentedLagrangian::setConstraintTolerance(double constraintTolerance) {
  epsilon = constraintTolerance;
}

double AugmentedLagrangian::getPenaltyStartValue() const { return mu0; }

void AugmentedLagrangian::setPenaltyStartValue(double penaltyStartValue) {
  mu0 = penaltyStartValue;
}

double AugmentedLagrangian::getPenaltyIncreaseFactor() const { return rhoMuPlus; }

void AugmentedLagrangian::setPenaltyIncreaseFactor(double penaltyIncreaseFactor) {
  rhoMuPlus = penaltyIncreaseFactor;
}

const base::DataMatrix& AugmentedLagrangian::getHistoryOfInnerIterationPoints() const {
  return xHistInner;
}

const std::vector<size_t>& AugmentedLagrangian::getHistoryOfInnerIterationNumbers() const {
  return kHistInner;
}

void AugmentedLagrangian::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new AugmentedLagrangian(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

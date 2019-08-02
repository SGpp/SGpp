// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/constrained/SquaredPenalty.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>

#include <vector>

namespace sgpp {
namespace optimization {
namespace optimizer {

namespace {
class PenalizedObjectiveFunction : public base::ScalarFunction {
 public:
  PenalizedObjectiveFunction(base::ScalarFunction& f, base::VectorFunction& g,
                             base::VectorFunction& h, double mu)
      : base::ScalarFunction(f.getNumberOfParameters()),
        f(f),
        g(g),
        h(h),
        mu(mu),
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
      if (gx[i] > 0.0) {
        value += mu * gx[i] * gx[i];
      }
    }

    for (size_t i = 0; i < mH; i++) {
      value += mu * hx[i] * hx[i];
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
  size_t mG;
  size_t mH;
};

class PenalizedObjectiveGradient : public base::ScalarFunctionGradient {
 public:
  PenalizedObjectiveGradient(base::ScalarFunctionGradient& fGradient,
                             base::VectorFunctionGradient& gGradient,
                             base::VectorFunctionGradient& hGradient, double mu)
      : base::ScalarFunctionGradient(fGradient.getNumberOfParameters()),
        fGradient(fGradient),
        gGradient(gGradient),
        hGradient(hGradient),
        mu(mu),
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

      if (gxi > 0.0) {
        value += mu * gxi * gxi;

        for (size_t t = 0; t < d; t++) {
          gradient[t] += mu * 2.0 * gxi * gradGx(i, t);
        }
      }
    }

    for (size_t i = 0; i < mH; i++) {
      const double hxi = hx[i];
      value += mu * hxi * hxi;

      for (size_t t = 0; t < d; t++) {
        gradient[t] += mu * 2.0 * hxi * gradHx(i, t);
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
  size_t mG;
  size_t mH;
};
}  // namespace

SquaredPenalty::SquaredPenalty(const base::ScalarFunction& f,
                               const base::VectorFunction& g,
                               const base::VectorFunction& h,
                               size_t maxItCount, double xTolerance, double constraintTolerance,
                               double penaltyStartValue, double penaltyIncreaseFactor)
    : ConstrainedOptimizer(f, g, h, maxItCount),
      theta(xTolerance),
      epsilon(constraintTolerance),
      mu0(penaltyStartValue),
      rhoMuPlus(penaltyIncreaseFactor),
      xHistInner(0, 0),
      kHistInner() {
}

SquaredPenalty::SquaredPenalty(const base::ScalarFunction& f,
                               const base::ScalarFunctionGradient& fGradient,
                               const base::VectorFunction& g,
                               const base::VectorFunctionGradient& gGradient,
                               const base::VectorFunction& h,
                               const base::VectorFunctionGradient& hGradient, size_t maxItCount,
                               double xTolerance, double constraintTolerance,
                               double penaltyStartValue, double penaltyIncreaseFactor)
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

SquaredPenalty::SquaredPenalty(const UnconstrainedOptimizer& unconstrainedOptimizer,
                               const base::VectorFunction& g,
                               const base::VectorFunctionGradient* gGradient,
                               const base::VectorFunction& h,
                               const base::VectorFunctionGradient* hGradient,
                               size_t maxItCount, double xTolerance, double constraintTolerance,
                               double penaltyStartValue, double penaltyIncreaseFactor)
    : ConstrainedOptimizer(unconstrainedOptimizer, g, gGradient, h, hGradient, maxItCount),
      theta(xTolerance),
      epsilon(constraintTolerance),
      mu0(penaltyStartValue),
      rhoMuPlus(penaltyIncreaseFactor),
      xHistInner(0, 0),
      kHistInner() {
}

SquaredPenalty::SquaredPenalty(const SquaredPenalty& other)
    : ConstrainedOptimizer(other),
      theta(other.theta),
      epsilon(other.epsilon),
      mu0(other.mu0),
      rhoMuPlus(other.rhoMuPlus),
      xHistInner(other.xHistInner),
      kHistInner(other.kHistInner) {
}

SquaredPenalty::~SquaredPenalty() {}

void SquaredPenalty::optimize() {
  base::Printer::getInstance().printStatusBegin("Optimizing (Squared Penalty)...");

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

  size_t breakIterationCounter = 0;
  const size_t BREAK_ITERATION_COUNTER_MAX = 10;
  size_t k = 1;

  const size_t unconstrainedN = N / 20;

  PenalizedObjectiveFunction fPenalized(*f, *g, *h, mu);
  std::unique_ptr<PenalizedObjectiveGradient> fPenalizedGradient;

  if ((fGradient != nullptr) && (gGradient != nullptr) && (hGradient != nullptr)) {
    fPenalizedGradient.reset(
        new PenalizedObjectiveGradient(*fGradient, *gGradient, *hGradient, mu));
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

double SquaredPenalty::getXTolerance() const { return theta; }

void SquaredPenalty::setXTolerance(double xTolerance) { theta = xTolerance; }

double SquaredPenalty::getConstraintTolerance() const { return epsilon; }

void SquaredPenalty::setConstraintTolerance(double constraintTolerance) {
  epsilon = constraintTolerance;
}

double SquaredPenalty::getPenaltyStartValue() const { return mu0; }

void SquaredPenalty::setPenaltyStartValue(double penaltyStartValue) { mu0 = penaltyStartValue; }

double SquaredPenalty::getPenaltyIncreaseFactor() const { return rhoMuPlus; }

void SquaredPenalty::setPenaltyIncreaseFactor(double penaltyIncreaseFactor) {
  rhoMuPlus = penaltyIncreaseFactor;
}

const base::DataMatrix& SquaredPenalty::getHistoryOfInnerIterationPoints() const {
  return xHistInner;
}

const std::vector<size_t>& SquaredPenalty::getHistoryOfInnerIterationNumbers() const {
  return kHistInner;
}

void SquaredPenalty::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new SquaredPenalty(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

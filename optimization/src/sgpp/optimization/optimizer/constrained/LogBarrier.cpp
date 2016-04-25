// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/constrained/LogBarrier.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>
#include <sgpp/optimization/function/vector/EmptyVectorFunction.hpp>

#include <vector>

namespace sgpp {
namespace optimization {
namespace optimizer {

namespace {
class PenalizedObjectiveFunction : public ScalarFunction {
 public:
  PenalizedObjectiveFunction(ScalarFunction& f, VectorFunction& g, double mu)
      : ScalarFunction(f.getNumberOfParameters()),
        f(f),
        g(g),
        mu(mu),
        m(g.getNumberOfComponents()) {}

  double eval(const base::DataVector& x) {
    for (size_t t = 0; t < d; t++) {
      if ((x[t] < 0.0) || (x[t] > 1.0)) {
        return INFINITY;
      }
    }

    const double fx = f.eval(x);
    base::DataVector gx(m);
    g.eval(x, gx);

    double value = fx;

    for (size_t i = 0; i < m; i++) {
      if (gx[i] < 0.0) {
        value -= mu * std::log(-gx[i]);
      } else {
        return INFINITY;
      }
    }

    return value;
  }

  void clone(std::unique_ptr<ScalarFunction>& clone) const {
    clone = std::unique_ptr<ScalarFunction>(new PenalizedObjectiveFunction(*this));
  }

  void setMu(double mu) { this->mu = mu; }

 protected:
  ScalarFunction& f;
  VectorFunction& g;
  double mu;
  size_t m;
};

class PenalizedObjectiveGradient : public ScalarFunctionGradient {
 public:
  PenalizedObjectiveGradient(ScalarFunctionGradient& fGradient, VectorFunctionGradient& gGradient,
                             double mu)
      : ScalarFunctionGradient(fGradient.getNumberOfParameters()),
        fGradient(fGradient),
        gGradient(gGradient),
        mu(mu),
        m(gGradient.getNumberOfComponents()) {}

  double eval(const base::DataVector& x, base::DataVector& gradient) {
    for (size_t t = 0; t < d; t++) {
      if ((x[t] < 0.0) || (x[t] > 1.0)) {
        gradient.setAll(NAN);
        return INFINITY;
      }
    }

    base::DataVector gradFx(d);
    const double fx = fGradient.eval(x, gradFx);
    base::DataVector gx(m);
    base::DataMatrix gradGx(m, d);
    gGradient.eval(x, gx, gradGx);

    double value = fx;
    gradient.resize(d);
    gradient = gradFx;

    for (size_t i = 0; i < m; i++) {
      const double gxi = gx[i];

      if (gxi < 0.0) {
        value -= mu * std::log(-gx[i]);

        for (size_t t = 0; t < d; t++) {
          gradient[t] -= mu * gradGx(i, t) / gx[i];
        }
      } else {
        return INFINITY;
      }
    }

    return value;
  }

  void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const {
    clone = std::unique_ptr<ScalarFunctionGradient>(new PenalizedObjectiveGradient(*this));
  }

  void setMu(double mu) { this->mu = mu; }

 protected:
  ScalarFunctionGradient& fGradient;
  VectorFunctionGradient& gGradient;
  double mu;
  size_t m;
};
}  // namespace

LogBarrier::LogBarrier(ScalarFunction& f, ScalarFunctionGradient& fGradient, VectorFunction& g,
                       VectorFunctionGradient& gGradient, size_t maxItCount, double tolerance,
                       double barrierStartValue, double barrierDecreaseFactor)
    : ConstrainedOptimizer(f, g, EmptyVectorFunction::getInstance(), maxItCount),
      fGradient(fGradient),
      gGradient(gGradient),
      theta(tolerance),
      mu0(barrierStartValue),
      rhoMuMinus(barrierDecreaseFactor),
      xHistInner(0, 0),
      kHistInner() {}

LogBarrier::~LogBarrier() {}

void LogBarrier::optimize() {
  Printer::getInstance().printStatusBegin("Optimizing (Log Barrier)...");

  const size_t d = f.getNumberOfParameters();

  xOpt.resize(0);
  fOpt = NAN;
  xHist.resize(0, d);
  fHist.resize(0);
  xHistInner.resize(0, d);
  kHistInner.clear();

  base::DataVector x(x0);
  double fx = f.eval(x);

  xHist.appendRow(x);
  fHist.append(fx);

  base::DataVector xNew(d);

  base::DataVector gx(g.getNumberOfComponents());

  double mu = mu0;

  size_t breakIterationCounter = 0;
  const size_t BREAK_ITERATION_COUNTER_MAX = 10;
  size_t k = 1;

  const size_t unconstrainedN = N / 20;

  PenalizedObjectiveFunction fPenalized(f, g, mu);
  PenalizedObjectiveGradient fPenalizedGradient(fGradient, gGradient, mu);

  while (k < N) {
    fPenalized.setMu(mu);
    fPenalizedGradient.setMu(mu);

    AdaptiveGradientDescent unconstrainedOptimizer(fPenalized, fPenalizedGradient, unconstrainedN,
                                                   10.0 * theta);
    unconstrainedOptimizer.setStartingPoint(x);
    unconstrainedOptimizer.optimize();
    xNew = unconstrainedOptimizer.getOptimalPoint();

    const base::DataMatrix& innerPoints = unconstrainedOptimizer.getHistoryOfOptimalPoints();
    const size_t numberInnerIterations = innerPoints.getNrows();
    k += numberInnerIterations;

    x = xNew;
    fx = f.eval(x);
    g.eval(x, gx);
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
    Printer::getInstance().printStatusUpdate(std::to_string(k) + " evaluations, x = " +
                                             x.toString() + ", f(x) = " + std::to_string(fx) +
                                             ", g(x) = " + gx.toString());

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
  fOpt = fx;
  Printer::getInstance().printStatusEnd();
}

ScalarFunctionGradient& LogBarrier::getObjectiveGradient() const { return fGradient; }

VectorFunctionGradient& LogBarrier::getInequalityConstraintGradient() const { return gGradient; }

double LogBarrier::getTolerance() const { return theta; }

void LogBarrier::setTolerance(double tolerance) { theta = tolerance; }

double LogBarrier::getBarrierStartValue() const { return mu0; }

void LogBarrier::setBarrierStartValue(double barrierStartValue) { mu0 = barrierStartValue; }

double LogBarrier::getBarrierDecreaseFactor() const { return rhoMuMinus; }

void LogBarrier::setBarrierDecreaseFactor(double barrierDecreaseFactor) {
  rhoMuMinus = barrierDecreaseFactor;
}

const base::DataMatrix& LogBarrier::getHistoryOfInnerIterationPoints() const {
  return xHistInner;
}

const std::vector<size_t>& LogBarrier::getHistoryOfInnerIterationNumbers() const {
  return kHistInner;
}

void LogBarrier::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new LogBarrier(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

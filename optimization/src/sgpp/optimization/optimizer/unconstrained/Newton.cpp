// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/system/FullSLE.hpp>
#include <sgpp/optimization/optimizer/unconstrained/LineSearchArmijo.hpp>
#include <sgpp/optimization/optimizer/unconstrained/Newton.hpp>

#include <algorithm>
#include <numeric>

namespace sgpp {
namespace optimization {
namespace optimizer {

Newton::Newton(const base::ScalarFunction& f, const base::ScalarFunctionHessian& fHessian,
               size_t max_it_count, double beta, double gamma, double tolerance, double epsilon,
               double alpha1, double alpha2, double p)
    : UnconstrainedOptimizer(f, nullptr, &fHessian, max_it_count),
      beta(beta),
      gamma(gamma),
      tol(tolerance),
      eps(epsilon),
      alpha1(alpha1),
      alpha2(alpha2),
      p(p),
      defaultSleSolver(base::sle_solver::GaussianElimination()),
      sleSolver(defaultSleSolver) {
}

Newton::Newton(const base::ScalarFunction& f, const base::ScalarFunctionHessian& fHessian,
               size_t max_it_count, double beta, double gamma, double tolerance, double epsilon,
               double alpha1, double alpha2, double p, const base::sle_solver::SLESolver& sleSolver)
    : UnconstrainedOptimizer(f, nullptr, &fHessian, max_it_count),
      beta(beta),
      gamma(gamma),
      tol(tolerance),
      eps(epsilon),
      alpha1(alpha1),
      alpha2(alpha2),
      p(p),
      defaultSleSolver(base::sle_solver::GaussianElimination()),
      sleSolver(sleSolver) {
}

Newton::Newton(const Newton& other)
    : UnconstrainedOptimizer(other),
      beta(other.beta),
      gamma(other.gamma),
      tol(other.tol),
      eps(other.eps),
      alpha1(other.alpha1),
      alpha2(other.alpha2),
      p(other.p),
      defaultSleSolver(base::sle_solver::GaussianElimination()),
      sleSolver(other.sleSolver) {
}

Newton::~Newton() {}

void Newton::optimize() {
  base::Printer::getInstance().printStatusBegin("Optimizing (Newton)...");

  const size_t d = f->getNumberOfParameters();

  xOpt.resize(0);
  fOpt = NAN;
  xHist.resize(0, d);
  fHist.resize(0);

  base::DataVector x(x0);
  double fx = NAN;

  base::DataVector gradFx(d);
  base::DataMatrix hessianFx(d, d);
  base::DataVector dk(d);
  base::DataVector s(d);
  base::DataVector y(d);

  FullSLE system(hessianFx);
  size_t k = 0;
  const bool statusPrintingEnabled = base::Printer::getInstance().isStatusPrintingEnabled();

  const double BINDING_TOLERANCE = 1e-6;

  while (k < N) {
    // calculate gradient, Hessian and gradient norm
    fx = fHessian->eval(x, gradFx, hessianFx);
    k++;

    if (k == 1) {
      xHist.appendRow(x);
      fHist.append(fx);
    }

    // RHS of linear system to be solved
    for (size_t t = 0; t < d; t++) {
      // search direction (negated gradient)
      s[t] = -gradFx[t];

      // is constraint binding?
      // (i.e., we are at the boundary and the search direction points outwards)
      if (((x[t] < BINDING_TOLERANCE) && (s[t] < 0.0)) ||
          ((x[t] > 1.0 - BINDING_TOLERANCE) && (s[t] > 0.0))) {
        // discard variable by setting RHS to zero
        s[t] = 0.0;
        gradFx[t] = 0.0;
        // eliminate variable from linear system matrix
        hessianFx(t, t) = 1.0;

        for (size_t t2 = 0; t2 < d; t2++) {
          if (t != t2) {
            hessianFx(t, t2) = 0.0;
            hessianFx(t2, t) = 0.0;
          }
        }
      }
    }

    const double sNorm = s.l2Norm();

    // exit if norm small enough
    if (sNorm < tol) {
      break;
    }

    // solve linear system with Hessian as system matrix
    if (statusPrintingEnabled) {
      base::Printer::getInstance().disableStatusPrinting();
    }

    bool lsSolved = sleSolver.solve(system, s, dk);

    if (statusPrintingEnabled) {
      base::Printer::getInstance().enableStatusPrinting();
    }

    // norm of solution
    const double dkNorm = dk.l2Norm();

    // acceptance criterion
    if (lsSolved &&
        (s.dotProduct(dk) >= std::min(alpha1, alpha2 * std::pow(dkNorm, p)) * dkNorm * dkNorm)) {
      // normalized solution as new search direction
      for (size_t t = 0; t < d; t++) {
        s[t] = dk[t] / dkNorm;
      }
    } else {
      // restart method (negated normalized gradient as new search direction)
      s.mult(1.0 / sNorm);
    }

    // status printing
    base::Printer::getInstance().printStatusUpdate(
        std::to_string(k) + " evaluations, x = " + x.toString() + ", f(x) = " + std::to_string(fx));

    // line search
    if (!lineSearchArmijo(*f, beta, gamma, tol, eps, x, fx, gradFx, s, y, k)) {
      // line search failed ==> exit
      // (either a "real" error occurred or the improvement
      // achieved is too small)
      break;
    }

    x = y;
    xHist.appendRow(x);
    fHist.append(fx);
  }

  xOpt.resize(d);
  xOpt = x;
  fOpt = fx;
  base::Printer::getInstance().printStatusEnd();
}

double Newton::getBeta() const { return beta; }

void Newton::setBeta(double beta) { this->beta = beta; }

double Newton::getGamma() const { return gamma; }

void Newton::setGamma(double gamma) { this->gamma = gamma; }

double Newton::getTolerance() const { return tol; }

void Newton::setTolerance(double tolerance) { tol = tolerance; }

double Newton::getEpsilon() const { return eps; }

void Newton::setEpsilon(double epsilon) { eps = epsilon; }

double Newton::getAlpha1() const { return alpha1; }

void Newton::setAlpha1(double alpha1) { this->alpha1 = alpha1; }

double Newton::getAlpha2() const { return alpha2; }

void Newton::setAlpha2(double alpha2) { this->alpha2 = alpha2; }

double Newton::getP() const { return p; }

void Newton::setP(double p) { this->p = p; }

void Newton::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new Newton(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

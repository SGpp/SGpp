// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/optimizer/unconstrained/GradientDescent.hpp>

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/LineSearchArmijo.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      GradientDescent::GradientDescent(
        ScalarFunction& f,
        ScalarFunctionGradient& fGradient,
        size_t maxItCount,
        float_t beta,
        float_t gamma,
        float_t tolerance,
        float_t epsilon) :
        UnconstrainedOptimizer(f, maxItCount),
        fGradient(fGradient),
        beta(beta),
        gamma(gamma),
        tol(tolerance),
        eps(epsilon) {
      }

      void GradientDescent::optimize() {
        printer.printStatusBegin("Optimizing (gradient descent)...");

        const size_t d = f.getDimension();

        xOpt.resize(0);
        fOpt = NAN;
        xHist.resize(0, d);
        fHist.resize(0);

        base::DataVector x(x0);
        float_t fx = NAN;

        base::DataVector gradFx(d);
        base::DataVector s(d);
        base::DataVector y(d);
        size_t k = 0;

        while (k < N) {
          // calculate gradient and norm
          fx = fGradient.eval(x, gradFx);
          k++;
          float_t gradFxNorm = gradFx.l2Norm();

          if (k == 1) {
            xHist.appendRow(x);
            fHist.append(fx);
          }

          // exit if norm small enough
          if (gradFxNorm < tol) {
            break;
          }

          // search direction is the normalized negated gradient
          for (size_t t = 0; t < d; t++) {
            s[t] = -gradFx[t] / gradFxNorm;
          }

          // status printing
          printer.printStatusUpdate(
            std::to_string(k) + " evaluations, x = " + x.toString() +
            ", f(x) = " + std::to_string(fx));

          // line search
          if (!lineSearchArmijo(f, beta, gamma, tol, eps, x, fx,
                                gradFx, s, y, k)) {
            // line search failed ==> exit
            // (either a "real" error occured or the improvement
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
        printer.printStatusEnd();
      }

      ScalarFunctionGradient& GradientDescent::getObjectiveGradient() const {
        return fGradient;
      }

      float_t GradientDescent::getBeta() const {
        return beta;
      }

      void GradientDescent::setBeta(float_t beta) {
        this->beta = beta;
      }

      float_t GradientDescent::getGamma() const {
        return gamma;
      }

      void GradientDescent::setGamma(float_t gamma) {
        this->gamma = gamma;
      }

      float_t GradientDescent::getTolerance() const {
        return tol;
      }

      void GradientDescent::setTolerance(float_t tolerance) {
        tol = tolerance;
      }

      float_t GradientDescent::getEpsilon() const {
        return eps;
      }

      void GradientDescent::setEpsilon(float_t epsilon) {
        eps = epsilon;
      }

      void GradientDescent::clone(
        std::unique_ptr<UnconstrainedOptimizer>& clone) const {
        clone = std::unique_ptr<UnconstrainedOptimizer>(
                  new GradientDescent(*this));
      }
    }
  }
}

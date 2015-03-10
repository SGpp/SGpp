// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/GradientMethod.hpp>
#include <sgpp/optimization/optimizer/LineSearchArmijo.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      GradientMethod::GradientMethod(
        ObjectiveFunction& f,
        ObjectiveGradient& fGradient,
        size_t N, float_t beta, float_t gamma, float_t tolerance,
        float_t epsilon) :
        Optimizer(f, N),
        fGradient(fGradient),
        beta(beta),
        gamma(gamma),
        tol(tolerance),
        eps(epsilon) {
      }

      float_t GradientMethod::optimize(std::vector<float_t>& xOpt) {
        printer.printStatusBegin("Optimizing (gradient method)...");

        const size_t d = f.getDimension();
        std::vector<float_t> x(x0);
        float_t fx = 0.0;

        base::DataVector gradFx(d);
        std::vector<float_t> s(d, 0.0);
        std::vector<float_t> y(d, 0.0);
        size_t k;

        for (k = 0; k < N; k++) {
          // calculate gradient and norm
          fx = fGradient.evalGradient(x, gradFx);
          float_t gradFxNorm = gradFx.l2Norm();

          // exit if norm small enough
          if (gradFxNorm < tol) {
            break;
          }

          // search direction is the normalized negated gradient
          for (size_t t = 0; t < d; t++) {
            s[t] = -gradFx[t] / gradFxNorm;
          }

          // status printing
          printer.printStatusUpdate(std::to_string(k) + " steps, f(x) = " +
                                    std::to_string(fx));

          // line search
          if (!lineSearchArmijo(f, beta, gamma, tol, eps, x, fx,
                                gradFx, s, y)) {
            // line search failed ==> exit
            // (either a "real" error occured or the improvement
            // achieved is too small)
            break;
          }

          x = y;
        }

        xOpt = x;

        printer.printStatusUpdate(std::to_string(k) + " steps, f(x) = " +
                                  std::to_string(fx));
        printer.printStatusEnd();

        return fx;
      }

      ObjectiveGradient& GradientMethod::getObjectiveGradient() const {
        return fGradient;
      }

      float_t GradientMethod::getBeta() const {
        return beta;
      }

      void GradientMethod::setBeta(float_t beta) {
        this->beta = beta;
      }

      float_t GradientMethod::getGamma() const {
        return gamma;
      }

      void GradientMethod::setGamma(float_t gamma) {
        this->gamma = gamma;
      }

      float_t GradientMethod::getTolerance() const {
        return tol;
      }

      void GradientMethod::setTolerance(float_t tolerance) {
        tol = tolerance;
      }

      float_t GradientMethod::getEpsilon() const {
        return eps;
      }

      void GradientMethod::setEpsilon(float_t epsilon) {
        eps = epsilon;
      }

    }
  }
}

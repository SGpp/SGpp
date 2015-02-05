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

      const float_t GradientMethod::DEFAULT_BETA = 0.5;
      const float_t GradientMethod::DEFAULT_GAMMA = 1e-2;
      const float_t GradientMethod::DEFAULT_TOLERANCE = 1e-8;
      const float_t GradientMethod::DEFAULT_EPSILON = 1e-18;

      GradientMethod::GradientMethod(
        const function::Objective& f,
        const function::ObjectiveGradient& fGradient,
        size_t N, float_t beta, float_t gamma, float_t tolerance, float_t epsilon) :
        Optimizer(f, N),
        fGradient(fGradient.clone()),
        beta(beta),
        gamma(gamma),
        tol(tolerance),
        eps(epsilon) {
      }

      float_t GradientMethod::optimize(std::vector<float_t>& xOpt) {
        tools::printer.printStatusBegin("Optimizing (gradient method)...");

        size_t d = f->getDimension();
        std::vector<float_t> x(x0);
        float_t fx = 0.0;

        base::DataVector gradFx(d);
        std::vector<float_t> s(d, 0.0);
        std::vector<float_t> y(d, 0.0);
        size_t k;

        for (k = 0; k < N; k++) {
          // calculate gradient and norm
          fx = fGradient->evalGradient(x, gradFx);
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
          tools::printer.printStatusUpdate(toString(k) + " steps, f(x) = " + toString(fx));

          // line search
          if (!lineSearchArmijo(*f, beta, gamma, tol, eps, x, fx, gradFx, s, y)) {
            // line search failed ==> exit
            // (either a "real" error occured or the improvement achieved is too small)
            break;
          }

          x = y;
        }

        xOpt = x;

        tools::printer.printStatusUpdate(toString(k) + " steps, f(x) = " + toString(fx));
        tools::printer.printStatusEnd();

        return fx;
      }

      Optimizer* GradientMethod::clone() {
        Optimizer* result =
          new GradientMethod(*f, *fGradient, N, beta, gamma, tol, eps);
        result->setStartingPoint(x0);
        return result;
      }

      function::ObjectiveGradient& GradientMethod::getObjectiveGradient() const {
        return *fGradient;
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

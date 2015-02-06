// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/NLCG.hpp>
#include <sgpp/optimization/optimizer/LineSearchArmijo.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <numeric>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      const float_t NLCG::DEFAULT_BETA = 0.5;
      const float_t NLCG::DEFAULT_GAMMA = 1e-2;
      const float_t NLCG::DEFAULT_TOLERANCE = 1e-8;
      const float_t NLCG::DEFAULT_EPSILON = 1e-18;
      const float_t NLCG::DEFAULT_RESTART_THRESHOLD = 0.1;

      NLCG::NLCG(const function::Objective& f, const function::ObjectiveGradient& fGradient,
                 size_t maxItCount, float_t beta, float_t gamma,
                 float_t tolerance, float_t epsilon, float_t restartThreshold) :
        Optimizer(f, maxItCount),
        beta(beta),
        gamma(gamma),
        tol(tolerance),
        eps(epsilon),
        alpha(restartThreshold) {
        function::ObjectiveGradient* fGradientPtr;
        fGradient.clone(fGradientPtr);
        this->fGradient.reset(fGradientPtr);
      }

      float_t NLCG::optimize(std::vector<float_t>& xOpt) {
        tools::printer.printStatusBegin("Optimizing (NLCG)...");

        size_t d = f->getDimension();
        std::vector<float_t> x(x0);
        float_t fx;
        float_t fy;

        base::DataVector gradFx(d);
        base::DataVector gradFy(d);
        std::vector<float_t> s(d, 0.0);
        std::vector<float_t> sNormalized(d, 0.0);
        std::vector<float_t> y(d, 0.0);
        size_t k;

        fx = fGradient->evalGradient(x0, gradFx);
        float_t gradFxNorm = gradFx.l2Norm();
        float_t gradFyNorm = 0.0;

        // negated gradient as starting search direction
        for (size_t t = 0; t < d; t++) {
          s[t] = -gradFx[t];
        }

        for (k = 0; k < N; k++) {
          // exit if norm small enough
          if (gradFxNorm < tol) {
            break;
          }

          // normalize search direction
          float_t s_norm = std::sqrt(std::inner_product(s.begin(), s.end(), s.begin(), 0.0));

          for (size_t t = 0; t < d; t++) {
            sNormalized[t] = s[t] / s_norm;
          }

          // line search
          if (!lineSearchArmijo(*f, beta, gamma, tol, eps, x, fx, gradFx, sNormalized, y)) {
            // line search failed ==> exit
            // (either a "real" error occured or the improvement achieved is too small)
            break;
          }

          // calculate gradient and norm
          fy = fGradient->evalGradient(y, gradFy);
          gradFyNorm = gradFy.l2Norm();

          float_t beta = 0.0;

          // the method is restarted (beta = 0), if the following criterion is *not* met
          if (std::abs(gradFy.dotProduct(gradFx)) / (gradFyNorm * gradFyNorm) < alpha) {
            // Polak-Ribiere coefficient
            for (size_t t = 0; t < d; t++) {
              beta += gradFy[t] * (gradFy[t] - gradFx[t]);
            }

            beta /= gradFxNorm * gradFxNorm;
          }

          // new search direction
          for (size_t t = 0; t < d; t++) {
            s[t] = beta * s[t] - gradFy[t];
          }

          // status printing
          tools::printer.printStatusUpdate(toString(k) + " steps, f(x) = " + toString(fx));

          x = y;
          fx = fy;
          gradFx = gradFy;
          gradFxNorm = gradFyNorm;
        }

        xOpt = x;

        tools::printer.printStatusUpdate(toString(k) + " steps, f(x) = " + toString(fx));
        tools::printer.printStatusEnd();

        return fx;
      }

      void NLCG::clone(Optimizer*& clone) {
        clone = new NLCG(*f, *fGradient, N, beta, gamma, tol, eps, alpha);
        clone->setStartingPoint(x0);
      }

      function::ObjectiveGradient& NLCG::getObjectiveGradient() const {
        return *fGradient;
      }

      float_t NLCG::getBeta() const {
        return beta;
      }

      void NLCG::setBeta(float_t beta) {
        this->beta = beta;
      }

      float_t NLCG::getGamma() const {
        return gamma;
      }

      void NLCG::setGamma(float_t gamma) {
        this->gamma = gamma;
      }

      float_t NLCG::getTolerance() const {
        return tol;
      }

      void NLCG::setTolerance(float_t tolerance) {
        tol = tolerance;
      }

      float_t NLCG::getEpsilon() const {
        return eps;
      }

      void NLCG::setEpsilon(float_t epsilon) {
        eps = epsilon;
      }

      float_t NLCG::getRestartThreshold() const {
        return alpha;
      }

      void NLCG::setRestartThreshold(float_t restartThreshold) {
        alpha = restartThreshold;
      }

    }
  }
}

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/NLCG.hpp>
#include <sgpp/optimization/optimizer/unconstrained/LineSearchArmijo.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <numeric>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      NLCG::NLCG(ObjectiveFunction& f,
                 ObjectiveGradient& fGradient,
                 size_t maxItCount, float_t beta, float_t gamma,
                 float_t tolerance, float_t epsilon,
                 float_t restartThreshold) :
        UnconstrainedOptimizer(f, maxItCount),
        fGradient(fGradient),
        beta(beta),
        gamma(gamma),
        tol(tolerance),
        eps(epsilon),
        alpha(restartThreshold) {
      }

      float_t NLCG::optimize(base::DataVector& xOpt) {
        printer.printStatusBegin("Optimizing (NLCG)...");

        const size_t d = f.getDimension();
        base::DataVector x(x0);
        float_t fx;
        float_t fy;

        base::DataVector gradFx(d);
        base::DataVector gradFy(d);
        base::DataVector s(d);
        base::DataVector sNormalized(d);
        base::DataVector y(d);

        fx = fGradient.eval(x0, gradFx);
        float_t gradFxNorm = gradFx.l2Norm();

        size_t k = 1;

        // negated gradient as starting search direction
        for (size_t t = 0; t < d; t++) {
          s[t] = -gradFx[t];
        }

        while (k < N) {
          // exit if norm small enough
          if (gradFxNorm < tol) {
            break;
          }

          // normalize search direction
          const float_t sNorm = s.l2Norm();

          for (size_t t = 0; t < d; t++) {
            sNormalized[t] = s[t] / sNorm;
          }

          // line search
          if (!lineSearchArmijo(f, beta, gamma, tol, eps, x, fx,
                                gradFx, sNormalized, y, k)) {
            // line search failed ==> exit
            // (either a "real" error occured or the improvement achieved is
            // too small)
            break;
          }

          // calculate gradient and norm
          fy = fGradient.eval(y, gradFy);
          k++;

          const float_t gradFyNorm = gradFy.l2Norm();

          float_t beta = 0.0;

          // the method is restarted (beta = 0), if the following criterion
          // is *not* met
          if (std::abs(gradFy.dotProduct(gradFx)) /
              (gradFyNorm * gradFyNorm) < alpha) {
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
          printer.printStatusUpdate(
            std::to_string(k) + " evaluations, f(x) = " +
            std::to_string(fx));

          x = y;
          fx = fy;
          gradFx = gradFy;
          gradFxNorm = gradFyNorm;
        }

        xOpt.resize(d);
        xOpt = x;

        printer.printStatusUpdate(
          std::to_string(k) + " evaluations, f(x) = " +
          std::to_string(fx));
        printer.printStatusEnd();

        return fx;
      }

      ObjectiveGradient& NLCG::getObjectiveGradient() const {
        return fGradient;
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

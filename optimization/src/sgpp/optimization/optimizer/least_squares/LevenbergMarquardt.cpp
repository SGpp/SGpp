// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/least_squares/LevenbergMarquardt.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      LevenbergMarquardt::LevenbergMarquardt(
        VectorFunction& phi,
        VectorFunctionGradient& phiGradient,
        size_t maxItCount,
        float_t tolerance,
        float_t initialDamping,
        float_t acceptanceThreshold,
        float_t effectivenessThreshold) :
        LeastSquaresOptimizer(phi, maxItCount),
        phiGradient(phiGradient),
        tol(tolerance),
        mu0(initialDamping),
        beta0(acceptanceThreshold),
        beta1(effectivenessThreshold),
        defaultSleSolver(sle_solver::GaussianElimination()),
        sleSolver(defaultSleSolver) {
      }

      LevenbergMarquardt::LevenbergMarquardt(
        VectorFunction& phi,
        VectorFunctionGradient& phiGradient,
        size_t maxItCount,
        float_t tolerance,
        float_t damping,
        float_t acceptanceThreshold,
        float_t effectivenessThreshold,
        const sle_solver::SLESolver& sleSolver) :
        LeastSquaresOptimizer(phi, maxItCount),
        phiGradient(phiGradient),
        tol(tolerance),
        mu0(damping),
        beta0(acceptanceThreshold),
        beta1(effectivenessThreshold),
        defaultSleSolver(sle_solver::GaussianElimination()),
        sleSolver(sleSolver) {
      }

      void LevenbergMarquardt::optimize() {
        Printer::getInstance().printStatusBegin("Optimizing (Levenberg-Marquardt)...");

        const size_t d = phi.getNumberOfParameters();
        const size_t m = phi.getNumberOfComponents();

        xOpt.resize(0);
        fOpt = NAN;
        xHist.resize(0, d);
        fHist.resize(0);

        base::DataVector x(x0);
        base::DataVector phix(m);
        base::DataMatrix gradPhix(m, d);

        base::DataVector y(d);
        base::DataVector phiy(m);

        base::DataMatrix A(d, d);
        FullSLE system(A);
        base::DataVector s(d);
        base::DataVector b(d);

        float_t fx = NAN;
        float_t mu = mu0;
        base::DataVector gradPhixTimesS(m);
        base::DataMatrix gradPhixSquared(m, m);

        size_t k = 0;
        const bool statusPrintingEnabled = Printer::getInstance().isStatusPrintingEnabled();

        while (k < N) {
          // calculate gradient of phi
          phiGradient.eval(x, phix, gradPhix);
          k++;

          // calculate f
          fx = 0.0;

          for (size_t i = 0; i < m; i++) {
            fx += std::pow(phix[i], 2.0);
          }

          if (k == 1) {
            xHist.appendRow(x);
            fHist.append(fx);
          }

          // calculate (nabla phi(x))' * (nabla phi(x))
          for (size_t t1 = 0; t1 < d; t1++) {
            for (size_t t2 = 0; t2 < d; t2++) {
              float_t entry = 0.0;

              for (size_t i = 0; i < m; i++) {
                entry += gradPhix(i, t1) * gradPhix(i, t2);
              }

              gradPhixSquared(t1, t2) = entry;
              A(t1, t2) = gradPhixSquared(t1, t2);
            }
          }

          // RHS of linear system to be solved
          for (size_t t = 0; t < d; t++) {
            float_t entry = 0.0;

            for (size_t i = 0; i < m; i++) {
              entry -= gradPhix(i, t) * phix[i];
            }

            b[t] = entry;
          }

          while (mu < 1e10) {
            // matrix of linear system to be solved
            for (size_t t = 0; t < d; t++) {
              A(t, t) = gradPhixSquared(t, t) + mu * mu;
            }

            // solve linear system
            if (statusPrintingEnabled) {
              Printer::getInstance().disableStatusPrinting();
            }

            bool lsSolved = sleSolver.solve(system, b, s);

            if (statusPrintingEnabled) {
              Printer::getInstance().enableStatusPrinting();
            }

            // fallback to mu * gradient of f, if linear system solving fails
            if (!lsSolved) {
              for (size_t t = 0; t < d; t++) {
                s[t] = -b[t];
              }

              const float_t sNorm = s.l2Norm();
              s.mult(mu / sNorm);
            }

            // calculate new point
            for (size_t t = 0; t < d; t++) {
              y[t] = x[t] + s[t];
            }

            if (!lsSolved) {
              mu *= 2.0;
              break;
            }

            // evaluate phi at new point
            phi.eval(y, phiy);
            k++;

            // evaluate f and linearized version of at new point
            float_t fy = 0.0;
            float_t fLinearY = 0.0;
            gradPhix.mult(s, gradPhixTimesS);

            for (size_t i = 0; i < m; i++) {
              fy += std::pow(phiy[i], 2.0);
              fLinearY += std::pow(phix[i] + gradPhixTimesS[i], 2.0);
            }

            // effectiveness parameter
            const float_t epsilon = (fx - fy) / (fx - fLinearY);

            if (epsilon <= beta0) {
              mu *= 2.0;
            } else if (epsilon < beta1) {
              break;
            } else {
              mu /= 2.0;
              break;
            }
          }

          // status printing
          Printer::getInstance().printStatusUpdate(
            std::to_string(k) + " evaluations, x = " + x.toString() +
            ", f(x) = " + std::to_string(fx));

          x = y;
          xHist.appendRow(x);
          fHist.append(fx);

          // stopping criterion
          if (s.l2Norm() < tol) {
            break;
          }
        }

        xOpt.resize(d);
        xOpt = x;
        fOpt = fx;
        Printer::getInstance().printStatusEnd();
      }

      VectorFunctionGradient& LevenbergMarquardt::getPhiGradient() const {
        return phiGradient;
      }

      float_t LevenbergMarquardt::getTolerance() const {
        return tol;
      }

      void LevenbergMarquardt::setTolerance(float_t tolerance) {
        tol = tolerance;
      }

      float_t LevenbergMarquardt::getInitialDamping() const {
        return mu0;
      }

      void LevenbergMarquardt::setInitialDamping(float_t initialDamping) {
        mu0 = initialDamping;
      }

      float_t LevenbergMarquardt::getAcceptanceThreshold() const {
        return beta0;
      }

      void LevenbergMarquardt::setAcceptanceThreshold(
        float_t acceptanceThreshold) {
        beta0 = acceptanceThreshold;
      }

      float_t LevenbergMarquardt::getEffectivenessThreshold() const {
        return beta1;
      }

      void LevenbergMarquardt::setEffectivenessThreshold(
        float_t effectivenessThreshold) {
        beta1 = effectivenessThreshold;
      }

      void LevenbergMarquardt::clone(
        std::unique_ptr<LeastSquaresOptimizer>& clone) const {
        clone = std::unique_ptr<LeastSquaresOptimizer>(
                  new LevenbergMarquardt(*this));
      }
    }
  }
}

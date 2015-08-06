// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/BiCGStab.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <cmath>
#include <numeric>

namespace SGPP {
  namespace optimization {
    namespace sle_solver {

      BiCGStab::BiCGStab() :
        BiCGStab(DEFAULT_MAX_IT_COUNT, DEFAULT_TOLERANCE,
                 base::DataVector(0)) {
      }

      BiCGStab::BiCGStab(size_t maxItCount, float_t tolerance,
                         const base::DataVector& x0) :
        SLESolver(),
        N(maxItCount),
        tol(tolerance),
        x0(x0) {
      }

      bool BiCGStab::solve(SLE& system, base::DataVector& b,
                           base::DataVector& x) const {
        printer.printStatusBegin("Solving linear system (BiCGStab)...");

        const size_t n = b.getSize();
        base::DataVector r(n, 0.0);

        if (n == 1) {
          const float_t A = system.getMatrixEntry(0, 0);

          if (A != 0.0) {
            x.resize(1);
            x[0] = b.get(0) / A;
            printer.printStatusEnd();
            return true;
          } else {
            printer.printStatusEnd("error: could not solve linear system!");
            return false;
          }
        }

        x.resize(n);

        if (x0.getSize() == n) {
          x = x0;
        } else {
          x.setAll(0.0);
        }

        system.matrixVectorMultiplication(x, r);

        for (size_t i = 0; i < n; i++) {
          r[i] = b.get(i) - r[i];
        }

        base::DataVector r0Hat(r);
        float_t rho = 1.0;
        float_t alpha = 1.0;
        float_t omega = 1.0;
        base::DataVector v(n, 0.0);
        base::DataVector p(n, 0.0);
        base::DataVector s(n, 0.0);
        base::DataVector t(n, 0.0);
        float_t rNormSquared = 0.0;
        size_t k = 0;

        for (k = 0; k < N; k++) {
          float_t last_rho = rho;
          rho = r0Hat.dotProduct(r);
          float_t beta = (rho / last_rho) * (alpha / omega);

          for (size_t i = 0; i < n; i++) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
          }

          system.matrixVectorMultiplication(p, v);
          alpha = rho / r0Hat.dotProduct(v);

          for (size_t i = 0; i < n; i++) {
            s[i] = r[i] - alpha * v[i];
          }

          system.matrixVectorMultiplication(s, t);
          omega = t.dotProduct(s) / t.dotProduct(t);

          if (std::isnan(omega)) {
            printer.printStatusEnd("error: could not solve linear system!");
            return false;
          }

          for (size_t i = 0; i < n; i++) {
            x[i] = x[i] + alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
          }

          rNormSquared = r.dotProduct(r);

          printer.printStatusUpdate("k = " + std::to_string(k) +
                                    ", residual norm = " +
                                    std::to_string(sqrt(rNormSquared)));

          if (rNormSquared < tol * tol) {
            break;
          }
        }

        printer.printStatusUpdate("k = " + std::to_string(k) +
                                  ", residual norm = " +
                                  std::to_string(sqrt(rNormSquared)));
        printer.printStatusEnd();
        return true;
      }

      size_t BiCGStab::getMaxItCount() const {
        return N;
      }

      void BiCGStab::setMaxItCount(size_t maxItCount) {
        N = maxItCount;
      }

      float_t BiCGStab::getTolerance() const {
        return tol;
      }

      void BiCGStab::setTolerance(float_t tolerance) {
        tol = tolerance;
      }

      const base::DataVector& BiCGStab::getStartingPoint() const {
        return x0;
      }

      void BiCGStab::setStartingPoint(const base::DataVector& startingPoint) {
        x0.resize(startingPoint.getSize());
        x0 = startingPoint;
      }

    }
  }
}

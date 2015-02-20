// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/Newton.hpp>
#include <sgpp/optimization/optimizer/LineSearchArmijo.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <algorithm>
#include <numeric>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      const float_t Newton::DEFAULT_BETA = 0.5;
      const float_t Newton::DEFAULT_GAMMA = 1e-2;
      const float_t Newton::DEFAULT_TOLERANCE = 1e-8;
      const float_t Newton::DEFAULT_EPSILON = 1e-18;
      const float_t Newton::DEFAULT_ALPHA1 = 1e-6;
      const float_t Newton::DEFAULT_ALPHA2 = 1e-6;
      const float_t Newton::DEFAULT_P = 0.1;

      Newton::Newton(
        const ObjectiveFunction& f,
        const ObjectiveHessian& fHessian,
        size_t max_it_count, float_t beta, float_t gamma,
        float_t tolerance, float_t epsilon, float_t alpha1,
        float_t alpha2, float_t p) :
        Optimizer(f, max_it_count),
        defaultSleSolver(sle_solver::BiCGStab()),
        sleSolver(defaultSleSolver) {
        initialize(fHessian, beta, gamma, tolerance, epsilon,
                   alpha1, alpha2, p);
      }

      Newton::Newton(
        const ObjectiveFunction& f,
        const ObjectiveHessian& fHessian,
        size_t max_it_count, float_t beta, float_t gamma,
        float_t tolerance, float_t epsilon, float_t alpha1,
        float_t alpha2, float_t p,
        const sle_solver::SLESolver& sleSolver) :
        Optimizer(f, max_it_count),
        defaultSleSolver(sle_solver::BiCGStab()),
        sleSolver(sleSolver) {
        initialize(fHessian, beta, gamma, tolerance, epsilon,
                   alpha1, alpha2, p);
      }

      void Newton::initialize(const ObjectiveHessian& fHessian,
                              float_t beta, float_t gamma,
                              float_t tolerance, float_t epsilon,
                              float_t alpha1, float_t alpha2, float_t p) {
        fHessian.clone(this->fHessian);
        this->beta = beta;
        this->gamma = gamma;
        this->tol = tolerance;
        this->eps = epsilon;
        this->alpha1 = alpha1;
        this->alpha2 = alpha2;
        this->p = p;
      }

      float_t Newton::optimize(std::vector<float_t>& xOpt) {
        printer.printStatusBegin("Optimizing (Newton)...");

        size_t d = f->getDimension();
        std::vector<float_t> x = x0;
        float_t fx = INFINITY;

        std::vector<float_t> lsSolverX0(d, 0.0);
        bool lsSolved;
        std::vector<float_t> dk(d, 0.0);

        base::DataVector grad_fx(d);
        base::DataMatrix hessianFx(d, d);
        std::vector<float_t> s(d, 0.0);
        std::vector<float_t> y(d, 0.0);

        FullSLE system(hessianFx);
        size_t k;

        for (k = 0; k < N; k++) {
          // calculate gradient, Hessian and gradient norm
          fx = fHessian->evalHessian(x, grad_fx, hessianFx);
          float_t gradFxNorm = grad_fx.l2Norm();

          // exit if norm small enough
          if (gradFxNorm < tol) {
            break;
          }

          // RHS of linear system to be solved
          for (size_t t = 0; t < d; t++) {
            s[t] = -grad_fx[t];
          }

          // solve linear system with Hessian as system matrix
          system.setA(hessianFx);
          printer.disableStatusPrinting();
          lsSolved = sleSolver.solve(system, s, dk);
          printer.enableStatusPrinting();

          // norm of solution
          float_t dkNorm = sqrt(std::inner_product(dk.begin(), dk.end(),
                                                   dk.begin(), 0.0));

          // acceptance criterion
          if (lsSolved && (std::inner_product(s.begin(), s.end(),
                                              dk.begin(), 0.0) >=
                           std::min(alpha1, alpha2 * std::pow(dkNorm, p)) *
                           dkNorm * dkNorm)) {
            // normalized solution as new search direction
            for (size_t t = 0; t < d; t++) {
              s[t] = dk[t] / dkNorm;
            }
          } else {
            // restart method
            // (negated normalized gradient as new search direction)
            for (size_t t = 0; t < d; t++) {
              s[t] = s[t] / gradFxNorm;
            }
          }

          // status printing
          printer.printStatusUpdate(std::to_string(k) + " steps, f(x) = " +
                                    std::to_string(fx));

          // line search
          if (!lineSearchArmijo(*f, beta, gamma, tol, eps, x, fx,
                                grad_fx, s, y)) {
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

      void Newton::clone(std::unique_ptr<Optimizer>& clone) const {
        clone = std::unique_ptr<Optimizer>(
                  new Newton(*f, *fHessian, N, beta, gamma, tol, eps,
                             alpha1, alpha2, p, sleSolver));
        clone->setStartingPoint(x0);
      }

      ObjectiveHessian& Newton::getObjectiveHessian() const {
        return *fHessian;
      }

      float_t Newton::getBeta() const {
        return beta;
      }

      void Newton::setBeta(float_t beta) {
        this->beta = beta;
      }

      float_t Newton::getGamma() const {
        return gamma;
      }

      void Newton::setGamma(float_t gamma) {
        this->gamma = gamma;
      }

      float_t Newton::getTolerance() const {
        return tol;
      }

      void Newton::setTolerance(float_t tolerance) {
        tol = tolerance;
      }

      float_t Newton::getEpsilon() const {
        return eps;
      }

      void Newton::setEpsilon(float_t epsilon) {
        eps = epsilon;
      }

      float_t Newton::getAlpha1() const {
        return alpha1;
      }

      void Newton::setAlpha1(float_t alpha1) {
        this->alpha1 = alpha1;
      }

      float_t Newton::getAlpha2() const {
        return alpha2;
      }

      void Newton::setAlpha2(float_t alpha2) {
        this->alpha2 = alpha2;
      }

      float_t Newton::getP() const {
        return p;
      }

      void Newton::setP(float_t p) {
        this->p = p;
      }

    }
  }
}

/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/optimizer/Newton.hpp"
#include "opt/optimizer/LineSearchArmijo.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "opt/sle/system/Full.hpp"
#include "opt/tools/Printer.hpp"

#include <algorithm>
#include <numeric>

namespace sg {
  namespace opt {
    namespace optimizer {

      const double Newton::DEFAULT_BETA = 0.5;
      const double Newton::DEFAULT_GAMMA = 1e-2;
      const double Newton::DEFAULT_TOLERANCE = 1e-8;
      const double Newton::DEFAULT_EPSILON = 1e-18;
      const double Newton::DEFAULT_ALPHA1 = 1e-6;
      const double Newton::DEFAULT_ALPHA2 = 1e-6;
      const double Newton::DEFAULT_P = 0.1;

      Newton::Newton(
        function::Objective& f,
        function::ObjectiveHessian& f_hessian,
        size_t max_it_count, double beta, double gamma,
        double tolerance, double epsilon, double alpha1, double alpha2, double p) :
        Optimizer(f, max_it_count),
        f_hessian(f_hessian.clone()),
        default_sle_solver(sle::solver::BiCGStab()),
        sle_solver(default_sle_solver) {
        initialize(beta, gamma, tolerance, epsilon, alpha1, alpha2, p);
      }

      Newton::Newton(
        function::Objective& f,
        function::ObjectiveHessian& f_hessian,
        size_t max_it_count, double beta, double gamma,
        double tolerance, double epsilon, double alpha1, double alpha2, double p,
        const sle::solver::Solver& sle_solver) :
        Optimizer(f, max_it_count),
        f_hessian(f_hessian.clone()),
        default_sle_solver(sle::solver::BiCGStab()),
        sle_solver(sle_solver) {
        initialize(beta, gamma, tolerance, epsilon, alpha1, alpha2, p);
      }

      void Newton::initialize(double beta, double gamma, double tolerance, double epsilon,
                              double alpha1, double alpha2, double p) {
        this->beta = beta;
        this->gamma = gamma;
        this->tol = tolerance;
        this->eps = epsilon;
        this->alpha1 = alpha1;
        this->alpha2 = alpha2;
        this->p = p;
      }

      double Newton::optimize(std::vector<double>& xopt) {
        tools::printer.printStatusBegin("Optimizing (Newton)...");

        size_t d = f->getDimension();
        std::vector<double> x = x0;
        double fx = INFINITY;

        std::vector<double> ls_solver_x0(d, 0.0);
        bool ls_solved;
        std::vector<double> dk(d, 0.0);

        base::DataVector grad_fx(d);
        base::DataMatrix hessian_fx(d, d);
        std::vector<double> s(d, 0.0);
        std::vector<double> y(d, 0.0);

        sle::system::Full system(hessian_fx);
        size_t k;

        for (k = 0; k < N; k++) {
          // calculate gradient, Hessian and gradient norm
          fx = f_hessian->evalHessian(x, grad_fx, hessian_fx);
          double grad_fx_norm = grad_fx.l2Norm();

          // exit if norm small enough
          if (grad_fx_norm < tol) {
            break;
          }

          // RHS of linear system to be solved
          for (size_t t = 0; t < d; t++) {
            s[t] = -grad_fx[t];
          }

          // solve linear system with Hessian as system matrix
          system.setA(hessian_fx);
          tools::printer.disableStatusPrinting();
          ls_solved = sle_solver.solve(system, s, dk);
          tools::printer.enableStatusPrinting();

          // norm of solution
          double dk_norm = sqrt(std::inner_product(dk.begin(), dk.end(), dk.begin(), 0.0));

          // acceptance criterion
          if (ls_solved && (std::inner_product(s.begin(), s.end(), dk.begin(), 0.0) >=
                            std::min(alpha1, alpha2*std::pow(dk_norm, p)) * dk_norm*dk_norm)) {
            // normalized solution as new search direction
            for (size_t t = 0; t < d; t++) {
              s[t] = dk[t] / dk_norm;
            }
          } else {
            // restart method (negated normalized gradient as new search direction)
            for (size_t t = 0; t < d; t++) {
              s[t] = s[t] / grad_fx_norm;
            }
          }

          // status printing
          tools::printer.printStatusUpdate(toString(k) + " steps, f(x) = " + toString(fx));

          // line search
          if (!lineSearchArmijo(*f, beta, gamma, tol, eps, x, fx, grad_fx, s, y)) {
            // line search failed ==> exit
            // (either a "real" error occured or the improvement achieved is too small)
            break;
          }

          x = y;
        }

        xopt = x;

        tools::printer.printStatusUpdate(toString(k) + " steps, f(x) = " + toString(fx));
        tools::printer.printStatusEnd();

        return fx;
      }

      tools::SmartPointer<Optimizer> Newton::clone() {
        tools::SmartPointer<Optimizer> result(new Newton(
                                                *f, *f_hessian, N, beta, gamma, tol, eps, alpha1, alpha2, p, sle_solver));
        result->setStartingPoint(x0);
        return result;
      }

      const tools::SmartPointer<function::ObjectiveHessian>& Newton::getObjectiveHessian() const {
        return f_hessian;
      }

      double Newton::getBeta() const {
        return beta;
      }

      void Newton::setBeta(double beta) {
        this->beta = beta;
      }

      double Newton::getGamma() const {
        return gamma;
      }

      void Newton::setGamma(double gamma) {
        this->gamma = gamma;
      }

      double Newton::getTolerance() const {
        return tol;
      }

      void Newton::setTolerance(double tolerance) {
        tol = tolerance;
      }

      double Newton::getEpsilon() const {
        return eps;
      }

      void Newton::setEpsilon(double epsilon) {
        eps = epsilon;
      }

      double Newton::getAlpha1() const {
        return alpha1;
      }

      void Newton::setAlpha1(double alpha1) {
        this->alpha1 = alpha1;
      }

      double Newton::getAlpha2() const {
        return alpha2;
      }

      void Newton::setAlpha2(double alpha2) {
        this->alpha2 = alpha2;
      }

      double Newton::getP() const {
        return p;
      }

      void Newton::setP(double p) {
        this->p = p;
      }

    }
  }
}

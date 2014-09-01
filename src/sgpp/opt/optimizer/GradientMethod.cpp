/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/optimizer/GradientMethod.hpp"
#include "opt/optimizer/LineSearchArmijo.hpp"
#include "opt/tools/Printer.hpp"

namespace sg {
  namespace opt {
    namespace optimizer {

      const double GradientMethod::DEFAULT_BETA = 0.5;
      const double GradientMethod::DEFAULT_GAMMA = 1e-2;
      const double GradientMethod::DEFAULT_TOLERANCE = 1e-8;
      const double GradientMethod::DEFAULT_EPSILON = 1e-18;

      GradientMethod::GradientMethod(
        function::Objective& f,
        function::ObjectiveGradient& f_gradient,
        size_t N, double beta, double gamma, double tolerance, double epsilon) :
        Optimizer(f, N),
        f_gradient(f_gradient.clone()),
        beta(beta),
        gamma(gamma),
        tol(tolerance),
        eps(epsilon) {
      }

      double GradientMethod::optimize(std::vector<double>& xopt) {
        tools::printer.printStatusBegin("Optimizing (gradient method)...");

        size_t d = f->getDimension();
        std::vector<double> x(x0);
        double fx = 0.0;

        base::DataVector grad_fx(d);
        std::vector<double> s(d, 0.0);
        std::vector<double> y(d, 0.0);
        size_t k;

        for (k = 0; k < N; k++) {
          // calculate gradient and norm
          fx = f_gradient->evalGradient(x, grad_fx);
          double grad_fx_norm = grad_fx.l2Norm();

          // exit if norm small enough
          if (grad_fx_norm < tol) {
            break;
          }

          // search direction is the normalized negated gradient
          for (size_t t = 0; t < d; t++) {
            s[t] = -grad_fx[t] / grad_fx_norm;
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

      tools::SmartPointer<Optimizer> GradientMethod::clone() {
        tools::SmartPointer<Optimizer> result(
          new GradientMethod(*f, *f_gradient, N, beta, gamma, tol, eps));
        result->setStartingPoint(x0);
        return result;
      }

      const tools::SmartPointer<function::ObjectiveGradient>
      & GradientMethod::getObjectiveGradient() const {
        return f_gradient;
      }

      double GradientMethod::getBeta() const {
        return beta;
      }

      void GradientMethod::setBeta(double beta) {
        this->beta = beta;
      }

      double GradientMethod::getGamma() const {
        return gamma;
      }

      void GradientMethod::setGamma(double gamma) {
        this->gamma = gamma;
      }

      double GradientMethod::getTolerance() const {
        return tol;
      }

      void GradientMethod::setTolerance(double tolerance) {
        tol = tolerance;
      }

      double GradientMethod::getEpsilon() const {
        return eps;
      }

      void GradientMethod::setEpsilon(double epsilon) {
        eps = epsilon;
      }

    }
  }
}

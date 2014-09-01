/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_OPTIMIZER_NEWTON_HPP
#define SGPP_OPT_OPTIMIZER_NEWTON_HPP

#include "opt/optimizer/Optimizer.hpp"
#include "opt/function/ObjectiveHessian.hpp"
#include "opt/sle/solver/Solver.hpp"
#include "opt/sle/solver/BiCGStab.hpp"

#include <cstddef>

namespace sg {
  namespace opt {
    namespace optimizer {

      /**
       * Gradient-based nonlinear conjugate gradient method.
       *
       * The method is restarted with the steepest descent direction
       * if the inner product of negated gradient and search direction is not big enough
       * (criterion depending on three parameters).
       */
      class Newton : public Optimizer {
        public:
          /// default beta (parameter for Armijo's rule)
          static const double DEFAULT_BETA;
          /// default gamma (parameter for Armijo's rule)
          static const double DEFAULT_GAMMA;
          /// default tolerance (parameter for Armijo's rule)
          static const double DEFAULT_TOLERANCE;
          /// default epsilon (parameter for Armijo's rule)
          static const double DEFAULT_EPSILON;
          /// default steepest descent restart parameter 1
          static const double DEFAULT_ALPHA1;
          /// default steepest descent restart parameter 2
          static const double DEFAULT_ALPHA2;
          /// default steepest descent restart exponent
          static const double DEFAULT_P;

          /**
           * Constructor.
           * By default, the BiCGStab method is used to solve the linear systems.
           *
           * @param f                 objective function
           * @param f_gradient        objective function gradient
           * @param max_it_count      maximal number of iterations
           * @param beta              beta (parameter for Armijo's rule)
           * @param gamma             gamma (parameter for Armijo's rule)
           * @param tolerance         tolerance (parameter for Armijo's rule)
           * @param epsilon           epsilon (parameter for Armijo's rule)
           * @param alpha1            steepest descent restart parameter 1
           * @param alpha2            steepest descent restart parameter 2
           * @param p                 steepest descent restart exponent
           */
          Newton(function::Objective& f,
                 function::ObjectiveHessian& f_hessian,
                 size_t max_it_count = DEFAULT_N,
                 double beta = DEFAULT_BETA,
                 double gamma = DEFAULT_GAMMA,
                 double tolerance = DEFAULT_TOLERANCE,
                 double epsilon = DEFAULT_EPSILON,
                 double alpha1 = DEFAULT_ALPHA1,
                 double alpha2 = DEFAULT_ALPHA2,
                 double p = DEFAULT_P);

          /**
           * Constructor.
           * Do not destruct the solver before this object!
           *
           * @param f                 objective function
           * @param f_gradient        objective function gradient
           * @param max_it_count      maximal number of iterations
           * @param beta              beta (parameter for Armijo's rule)
           * @param gamma             gamma (parameter for Armijo's rule)
           * @param tolerance         tolerance (parameter for Armijo's rule)
           * @param epsilon           epsilon (parameter for Armijo's rule)
           * @param alpha1            steepest descent restart parameter 1
           * @param alpha2            steepest descent restart parameter 2
           * @param p                 steepest descent restart exponent
           * @param sle_solver        reference to linear solver for solving the linear systems
           *                          (Hessian as coefficient matrix)
           */
          Newton(function::Objective& f,
                 function::ObjectiveHessian& f_hessian,
                 size_t max_it_count,
                 double beta,
                 double gamma,
                 double tolerance,
                 double epsilon,
                 double alpha1,
                 double alpha2,
                 double p,
                 const sle::solver::Solver& sle_solver);

          /**
           * @param[out] xopt optimal point
           * @return          optimal objective function value
           */
          double optimize(std::vector<double>& xopt);

          /**
           * @return smart pointer to cloned object
           */
          tools::SmartPointer<Optimizer> clone();

          /**
           * @return smart pointer to objective function Hessian
           */
          const tools::SmartPointer<function::ObjectiveHessian>& getObjectiveHessian() const;

          /**
           * @return              beta (parameter for Armijo's rule)
           */
          double getBeta() const;

          /**
           * @param beta          beta (parameter for Armijo's rule)
           */
          void setBeta(double beta);

          /**
           * @return              gamma (parameter for Armijo's rule)
           */
          double getGamma() const;

          /**
           * @param gamma         gamma (parameter for Armijo's rule)
           */
          void setGamma(double gamma);

          /**
           * @return              tolerance (parameter for Armijo's rule)
           */
          double getTolerance() const;

          /**
           * @param tolerance     tolerance (parameter for Armijo's rule)
           */
          void setTolerance(double tolerance);

          /**
           * @return              epsilon (parameter for Armijo's rule)
           */
          double getEpsilon() const;

          /**
           * @param epsilon       epsilon (parameter for Armijo's rule)
           */
          void setEpsilon(double epsilon);

          /**
           * @return          steepest descent restart parameter 1
           */
          double getAlpha1() const;

          /**
           * @param alpha1    steepest descent restart parameter 1
           */
          void setAlpha1(double alpha1);

          /**
           * @return          steepest descent restart parameter 2
           */
          double getAlpha2() const;

          /**
           * @param alpha2    steepest descent restart parameter 2
           */
          void setAlpha2(double alpha2);

          /**
           * @return          steepest descent restart exponent
           */
          double getP() const;

          /**
           * @param p         steepest descent restart exponent
           */
          void setP(double p);

        protected:
          /// objective function Hessian
          tools::SmartPointer<function::ObjectiveHessian> f_hessian;
          /// beta (parameter for Armijo's rule)
          double beta;
          /// gamma (parameter for Armijo's rule)
          double gamma;
          /// tolerance (parameter for Armijo's rule)
          double tol;
          /// epsilon (parameter for Armijo's rule)
          double eps;
          /// steepest descent restart parameter 1
          double alpha1;
          /// steepest descent restart parameter 2
          double alpha2;
          /// steepest descent restart exponent
          double p;
          /// default linear solver
          const sle::solver::BiCGStab default_sle_solver;
          /// linear solver
          const sle::solver::Solver& sle_solver;

          /**
           * Internal function for initializing the member variables.
           *
           * @param beta              beta (parameter for Armijo's rule)
           * @param gamma             gamma (parameter for Armijo's rule)
           * @param tolerance         tolerance (parameter for Armijo's rule)
           * @param epsilon           epsilon (parameter for Armijo's rule)
           * @param alpha1            steepest descent restart parameter 1
           * @param alpha2            steepest descent restart parameter 2
           * @param p                 steepest descent restart exponent
           */
          void initialize(double beta, double gamma, double tolerance, double epsilon,
                          double alpha1, double alpha2, double p);
      };

    }
  }
}

#endif

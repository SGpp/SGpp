// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_NEWTON_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_NEWTON_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/Optimizer.hpp>
#include <sgpp/optimization/function/ObjectiveHessian.hpp>
#include <sgpp/optimization/sle/solver/BiCGStab.hpp>
#include <sgpp/optimization/sle/solver/SLESolver.hpp>

#include <cstddef>
#include <memory>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Gradient-based nonlinear conjugate gradient method.
       *
       * The method is restarted with the steepest descent direction
       * if the inner product of negated gradient and search direction is
       * not big enough (criterion depending on three parameters).
       */
      class Newton : public Optimizer {
        public:
          /// default beta (parameter for Armijo's rule)
          static constexpr float_t DEFAULT_BETA = 0.5;
          /// default gamma (parameter for Armijo's rule)
          static constexpr float_t DEFAULT_GAMMA = 1e-2;
          /// default tolerance (parameter for Armijo's rule)
          static constexpr float_t DEFAULT_TOLERANCE = 1e-8;
          /// default epsilon (parameter for Armijo's rule)
          static constexpr float_t DEFAULT_EPSILON = 1e-18;
          /// default steepest descent restart parameter 1
          static constexpr float_t DEFAULT_ALPHA1 = 1e-6;
          /// default steepest descent restart parameter 2
          static constexpr float_t DEFAULT_ALPHA2 = 1e-6;
          /// default steepest descent restart exponent
          static constexpr float_t DEFAULT_P = 0.1;

          /**
           * Constructor.
           * By default, the BiCGStab method is used to solve the
           * linear systems.
           *
           * @param f                 objective function
           * @param fHessian          objective function Hessian
           * @param maxItCount        maximal number of iterations
           * @param beta              beta (parameter for Armijo's rule)
           * @param gamma             gamma (parameter for Armijo's rule)
           * @param tolerance         tolerance (parameter for Armijo's rule)
           * @param epsilon           epsilon (parameter for Armijo's rule)
           * @param alpha1            steepest descent restart parameter 1
           * @param alpha2            steepest descent restart parameter 2
           * @param p                 steepest descent restart exponent
           */
          Newton(ObjectiveFunction& f,
                 ObjectiveHessian& fHessian,
                 size_t maxItCount = DEFAULT_N,
                 float_t beta = DEFAULT_BETA,
                 float_t gamma = DEFAULT_GAMMA,
                 float_t tolerance = DEFAULT_TOLERANCE,
                 float_t epsilon = DEFAULT_EPSILON,
                 float_t alpha1 = DEFAULT_ALPHA1,
                 float_t alpha2 = DEFAULT_ALPHA2,
                 float_t p = DEFAULT_P);

          /**
           * Constructor.
           * Do not destruct the solver before this object!
           *
           * @param f                 objective function
           * @param fHessian          objective function Hessian
           * @param maxItCount        maximal number of iterations
           * @param beta              beta (parameter for Armijo's rule)
           * @param gamma             gamma (parameter for Armijo's rule)
           * @param tolerance         tolerance (parameter for Armijo's rule)
           * @param epsilon           epsilon (parameter for Armijo's rule)
           * @param alpha1            steepest descent restart parameter 1
           * @param alpha2            steepest descent restart parameter 2
           * @param p                 steepest descent restart exponent
           * @param sleSolver         reference to linear solver for solving
           *                          the linear systems
           *                          (Hessian as coefficient matrix)
           */
          Newton(ObjectiveFunction& f,
                 ObjectiveHessian& fHessian,
                 size_t maxItCount,
                 float_t beta,
                 float_t gamma,
                 float_t tolerance,
                 float_t epsilon,
                 float_t alpha1,
                 float_t alpha2,
                 float_t p,
                 const sle_solver::SLESolver& sleSolver);

          /**
           * @param[out] xOpt optimal point
           * @return          optimal objective function value
           */
          float_t optimize(std::vector<float_t>& xOpt);

          /**
           * @return objective function Hessian
           */
          ObjectiveHessian& getObjectiveHessian() const;

          /**
           * @return              beta (parameter for Armijo's rule)
           */
          float_t getBeta() const;

          /**
           * @param beta          beta (parameter for Armijo's rule)
           */
          void setBeta(float_t beta);

          /**
           * @return              gamma (parameter for Armijo's rule)
           */
          float_t getGamma() const;

          /**
           * @param gamma         gamma (parameter for Armijo's rule)
           */
          void setGamma(float_t gamma);

          /**
           * @return              tolerance (parameter for Armijo's rule)
           */
          float_t getTolerance() const;

          /**
           * @param tolerance     tolerance (parameter for Armijo's rule)
           */
          void setTolerance(float_t tolerance);

          /**
           * @return              epsilon (parameter for Armijo's rule)
           */
          float_t getEpsilon() const;

          /**
           * @param epsilon       epsilon (parameter for Armijo's rule)
           */
          void setEpsilon(float_t epsilon);

          /**
           * @return          steepest descent restart parameter 1
           */
          float_t getAlpha1() const;

          /**
           * @param alpha1    steepest descent restart parameter 1
           */
          void setAlpha1(float_t alpha1);

          /**
           * @return          steepest descent restart parameter 2
           */
          float_t getAlpha2() const;

          /**
           * @param alpha2    steepest descent restart parameter 2
           */
          void setAlpha2(float_t alpha2);

          /**
           * @return          steepest descent restart exponent
           */
          float_t getP() const;

          /**
           * @param p         steepest descent restart exponent
           */
          void setP(float_t p);

        protected:
          /// objective function Hessian
          ObjectiveHessian& fHessian;
          /// beta (parameter for Armijo's rule)
          float_t beta;
          /// gamma (parameter for Armijo's rule)
          float_t gamma;
          /// tolerance (parameter for Armijo's rule)
          float_t tol;
          /// epsilon (parameter for Armijo's rule)
          float_t eps;
          /// steepest descent restart parameter 1
          float_t alpha1;
          /// steepest descent restart parameter 2
          float_t alpha2;
          /// steepest descent restart exponent
          float_t p;
          /// default linear solver
          const sle_solver::BiCGStab defaultSleSolver;
          /// linear solver
          const sle_solver::SLESolver& sleSolver;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_NEWTON_HPP */

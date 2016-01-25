// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NLCG_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NLCG_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Gradient-based nonlinear conjugate gradient method.
       *
       * The method is restarted with the steepest descent direction
       * if the normalized absolute value of the inner product of
       * two successive gradients
       * exceeds a "restart threshold" \f$\alpha\f$.
       */
      class NLCG : public UnconstrainedOptimizer {
        public:
          /// default beta (parameter for Armijo's rule)
          static constexpr float_t DEFAULT_BETA = 0.5;
          /// default gamma (parameter for Armijo's rule)
          static constexpr float_t DEFAULT_GAMMA = 1e-2;
          /// default tolerance (parameter for Armijo's rule)
          static constexpr float_t DEFAULT_TOLERANCE = 1e-8;
          /// default epsilon (parameter for Armijo's rule)
          static constexpr float_t DEFAULT_EPSILON = 1e-18;
          /// default restart threshold
          static constexpr float_t DEFAULT_RESTART_THRESHOLD = 0.1;

          /**
           * Constructor.
           *
           * @param f                 objective function
           * @param fGradient         objective function gradient
           * @param maxItCount        maximal number of iterations
           * @param beta              beta (parameter for Armijo's rule)
           * @param gamma             gamma (parameter for Armijo's rule)
           * @param tolerance         tolerance (parameter for Armijo's rule)
           * @param epsilon           epsilon (parameter for Armijo's rule)
           * @param restartThreshold  restart threshold
           */
          NLCG(ScalarFunction& f,
               ScalarFunctionGradient& fGradient,
               size_t maxItCount = DEFAULT_N,
               float_t beta = DEFAULT_BETA,
               float_t gamma = DEFAULT_GAMMA,
               float_t tolerance = DEFAULT_TOLERANCE,
               float_t epsilon = DEFAULT_EPSILON,
               float_t restartThreshold = DEFAULT_RESTART_THRESHOLD);

          /**
           * Destructor.
           */
          virtual ~NLCG() override;

          virtual void optimize() override;

          /**
           * @return objective function gradient
           */
          ScalarFunctionGradient& getObjectiveGradient() const;

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
           * @return                  restart threshold
           */
          float_t getRestartThreshold() const;

          /**
           * @param restartThreshold  restart threshold
           */
          void setRestartThreshold(float_t restartThreshold);

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const override;

        protected:
          /// objective function gradient
          ScalarFunctionGradient& fGradient;
          /// beta (parameter for Armijo's rule)
          float_t beta;
          /// gamma (parameter for Armijo's rule)
          float_t gamma;
          /// tolerance (parameter for Armijo's rule)
          float_t tol;
          /// epsilon (parameter for Armijo's rule)
          float_t eps;
          /// restart threshold
          float_t alpha;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NLCG_HPP */

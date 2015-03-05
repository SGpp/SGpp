// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_NLCG_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_NLCG_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/Optimizer.hpp>
#include <sgpp/optimization/function/ObjectiveGradient.hpp>

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
      class NLCG : public Optimizer {
        public:
          /// default beta (parameter for Armijo's rule)
          static const float_t DEFAULT_BETA;
          /// default gamma (parameter for Armijo's rule)
          static const float_t DEFAULT_GAMMA;
          /// default tolerance (parameter for Armijo's rule)
          static const float_t DEFAULT_TOLERANCE;
          /// default epsilon (parameter for Armijo's rule)
          static const float_t DEFAULT_EPSILON;
          /// default restart threshold
          static const float_t DEFAULT_RESTART_THRESHOLD;

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
          NLCG(const ObjectiveFunction& f,
               const ObjectiveGradient& fGradient,
               size_t maxItCount = DEFAULT_N,
               float_t beta = DEFAULT_BETA,
               float_t gamma = DEFAULT_GAMMA,
               float_t tolerance = DEFAULT_TOLERANCE,
               float_t epsilon = DEFAULT_EPSILON,
               float_t restartThreshold = DEFAULT_RESTART_THRESHOLD);

          /**
           * @param[out] xOpt optimal point
           * @return          optimal objective function value
           */
          float_t optimize(std::vector<float_t>& xOpt);

          /**
           * @param[out] clone pointer to cloned object
           */
          void clone(std::unique_ptr<Optimizer>& clone) const;

          /**
           * @return objective function gradient
           */
          ObjectiveGradient& getObjectiveGradient() const;

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

        protected:
          /// objective function gradient
          std::unique_ptr<ObjectiveGradient> fGradient;
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

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_NLCG_HPP */

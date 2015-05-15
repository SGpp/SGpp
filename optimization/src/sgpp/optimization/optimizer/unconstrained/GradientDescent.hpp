// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_GRADIENTDESCENT_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_GRADIENTDESCENT_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/ObjectiveGradient.hpp>

#include <memory>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Gradient-based method of steepest descent.
       */
      class GradientDescent : public UnconstrainedOptimizer {
        public:
          /// default maximal number of iterations
          static const size_t DEFAULT_MAX_IT_COUNT = 2000;
          /// default beta (parameter for Armijo's rule)
          static constexpr float_t DEFAULT_BETA = 0.5;
          /// default gamma (parameter for Armijo's rule)
          static constexpr float_t DEFAULT_GAMMA = 1e-2;
          /// default tolerance (parameter for Armijo's rule)
          static constexpr float_t DEFAULT_TOLERANCE = 1e-8;
          /// default epsilon (parameter for Armijo's rule)
          static constexpr float_t DEFAULT_EPSILON = 1e-18;

          /**
           * Constructor.
           *
           * @param f             objective function
           * @param fGradient     objective function gradient
           * @param maxItCount    maximal number of iterations
           * @param beta          beta (parameter for Armijo's rule)
           * @param gamma         gamma (parameter for Armijo's rule)
           * @param tolerance     tolerance (parameter for Armijo's rule)
           * @param epsilon       epsilon (parameter for Armijo's rule)
           */
          GradientDescent(ObjectiveFunction& f,
                          ObjectiveGradient& fGradient,
                          size_t maxItCount = DEFAULT_MAX_IT_COUNT,
                          float_t beta = DEFAULT_BETA,
                          float_t gamma = DEFAULT_GAMMA,
                          float_t tolerance = DEFAULT_TOLERANCE,
                          float_t epsilon = DEFAULT_EPSILON);

          /**
           * @param[out] xOpt optimal point
           * @return          optimal objective function value
           */
          float_t optimize(base::DataVector& xOpt);

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

        protected:
          /// objective function gradient
          ObjectiveGradient& fGradient;
          /// beta (parameter for Armijo's rule)
          float_t beta;
          /// gamma (parameter for Armijo's rule)
          float_t gamma;
          /// tolerance (parameter for Armijo's rule)
          float_t tol;
          /// epsilon (parameter for Armijo's rule)
          float_t eps;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_GRADIENTDESCENT_HPP */

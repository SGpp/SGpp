// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_RPROP_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_RPROP_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Rprop method for unconstrained optimization.
       */
      class Rprop : public UnconstrainedOptimizer {
        public:
          /// default tolerance
          static constexpr float_t DEFAULT_TOLERANCE = 1e-6;
          /// default initial step size
          static constexpr float_t DEFAULT_INITIAL_STEP_SIZE = 0.01;
          /// default step size increase factor
          static constexpr float_t DEFAULT_STEP_SIZE_INCREASE_FACTOR = 1.2;
          /// default step size decrease factor
          static constexpr float_t DEFAULT_STEP_SIZE_DECREASE_FACTOR = 0.5;

          /**
           * Constructor.
           *
           * @param f                       objective function
           * @param fGradient               objective function gradient
           * @param maxItCount              maximal number of
           *                                function evaluations
           * @param tolerance               tolerance
           * @param initialStepSize         initial step size
           * @param stepSizeIncreaseFactor  step size increase factor
           * @param stepSizeDecreaseFactor  step size decrease factor
           */
          Rprop(ObjectiveFunction& f,
                ObjectiveGradient& fGradient,
                size_t maxItCount = DEFAULT_N,
                float_t tolerance = DEFAULT_TOLERANCE,
                float_t initialStepSize = DEFAULT_INITIAL_STEP_SIZE,
                float_t stepSizeIncreaseFactor =
                  DEFAULT_STEP_SIZE_INCREASE_FACTOR,
                float_t stepSizeDecreaseFactor =
                  DEFAULT_STEP_SIZE_DECREASE_FACTOR);

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
           * @return tolerance
           */
          float_t getTolerance() const;

          /**
           * @param tolerance tolerance
           */
          void setTolerance(float_t tolerance);

          /**
           * @return initial step size
           */
          float_t getInitialStepSize() const;

          /**
           * @param initialStepSize initial step size
           */
          void setInitialStepSize(float_t initialStepSize);

          /**
           * @return step size increase factor
           */
          float_t getStepSizeIncreaseFactor() const;

          /**
           * @param stepSizeIncreaseFactor step size increase factor
           */
          void setStepSizeIncreaseFactor(float_t stepSizeIncreaseFactor);

          /**
           * @return step size decrease factor
           */
          float_t getStepSizeDecreaseFactor() const;

          /**
           * @param stepSizeDecreaseFactor step size decrease factor
           */
          void setStepSizeDecreaseFactor(float_t stepSizeDecreaseFactor);

        protected:
          /// objective function gradient
          ObjectiveGradient& fGradient;
          /// tolerance
          float_t theta;
          /// initial step size
          float_t initialAlpha;
          /// step size increase factor
          float_t rhoAlphaPlus;
          /// step size decrease factor
          float_t rhoAlphaMinus;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_RPROP_HPP */

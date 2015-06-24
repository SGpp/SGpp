// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_ADAPTIVEGRADIENTDESCENT_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_ADAPTIVEGRADIENTDESCENT_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/ObjectiveGradient.hpp>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Gradient descent with adaptive step size.
       */
      class AdaptiveGradientDescent : public UnconstrainedOptimizer {
        public:
          /// default tolerance
          static constexpr float_t DEFAULT_TOLERANCE = 1e-6;
          /// default step size increase factor
          static constexpr float_t DEFAULT_STEP_SIZE_INCREASE_FACTOR = 1.2;
          /// default step size decrease factor
          static constexpr float_t DEFAULT_STEP_SIZE_DECREASE_FACTOR = 0.5;
          /// default line search accuracy
          static constexpr float_t DEFAULT_LINE_SEARCH_ACCURACY = 0.01;

          /**
           * Constructor.
           *
           * @param f                       objective function
           * @param fGradient               objective function gradient
           * @param maxItCount              maximal number of
           *                                function evaluations
           * @param tolerance               tolerance
           * @param stepSizeIncreaseFactor  step size increase factor
           * @param stepSizeDecreaseFactor  step size decrease factor
           * @param lineSearchAccuracy      line search accuracy
           */
          AdaptiveGradientDescent(ObjectiveFunction& f,
                                  ObjectiveGradient& fGradient,
                                  size_t maxItCount = DEFAULT_N,
                                  float_t tolerance = DEFAULT_TOLERANCE,
                                  float_t stepSizeIncreaseFactor =
                                    DEFAULT_STEP_SIZE_INCREASE_FACTOR,
                                  float_t stepSizeDecreaseFactor =
                                    DEFAULT_STEP_SIZE_DECREASE_FACTOR,
                                  float_t lineSearchAccuracy = DEFAULT_LINE_SEARCH_ACCURACY);

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

          /**
           * @return line search accuracy
           */
          float_t getLineSearchAccuracy() const;

          /**
           * @param lineSearchAccuracy line search accuracy
           */
          void setLineSearchAccuracy(float_t lineSearchAccuracy);

        protected:
          /// objective function gradient
          ObjectiveGradient& fGradient;
          /// tolerance
          float_t theta;
          /// step size increase factor
          float_t rhoAlphaPlus;
          /// step size decrease factor
          float_t rhoAlphaMinus;
          /// line search accuracy
          float_t rhoLs;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_ADAPTIVEGRADIENTDESCENT_HPP */

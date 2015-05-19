// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_LOGBARRIER_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_LOGBARRIER_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/ObjectiveGradient.hpp>
#include <sgpp/optimization/function/ConstraintGradient.hpp>
#include <sgpp/optimization/optimizer/constrained/ConstrainedOptimizer.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      class LogBarrier : public ConstrainedOptimizer {
        public:
          static constexpr float_t DEFAULT_TOLERANCE = 1e-6;
          static constexpr float_t DEFAULT_BARRIER_START_VALUE = 1.0;
          static constexpr float_t DEFAULT_BARRIER_DECREASE_FACTOR = 0.5;

          LogBarrier(ObjectiveFunction& f,
                     ObjectiveGradient& fGradient,
                     ConstraintFunction& g,
                     ConstraintGradient& gGradient,
                     size_t maxItCount = DEFAULT_N,
                     float_t tolerance = DEFAULT_TOLERANCE,
                     float_t barrierStartValue = DEFAULT_BARRIER_START_VALUE,
                     float_t barrierDecreaseFactor =
                       DEFAULT_BARRIER_DECREASE_FACTOR);

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
           * @return inequality constraint function gradient
           */
          ConstraintGradient& getInequalityConstraintGradient() const;

          float_t getTolerance() const;
          void setTolerance(float_t tolerance);

          float_t getBarrierStartValue() const;
          void setBarrierStartValue(float_t barrierStartValue);

          float_t getBarrierDecreaseFactor() const;
          void setBarrierDecreaseFactor(float_t barrierDecreaseFactor);

        protected:
          /// objective function gradient
          ObjectiveGradient& fGradient;
          /// inequality constraint function gradient
          ConstraintGradient& gGradient;
          float_t theta;
          float_t mu0;
          float_t rhoMuMinus;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_LOGBARRIER_HPP */

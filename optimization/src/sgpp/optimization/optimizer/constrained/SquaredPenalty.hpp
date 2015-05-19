// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_SQUAREDPENALTY_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_SQUAREDPENALTY_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/ObjectiveGradient.hpp>
#include <sgpp/optimization/function/ConstraintGradient.hpp>
#include <sgpp/optimization/optimizer/constrained/ConstrainedOptimizer.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      class SquaredPenalty : public ConstrainedOptimizer {
        public:
          static constexpr float_t DEFAULT_X_TOLERANCE = 1e-6;
          static constexpr float_t DEFAULT_CONSTRAINT_TOLERANCE = 1e-6;
          static constexpr float_t DEFAULT_PENALTY_START_VALUE = 1.0;
          static constexpr float_t DEFAULT_PENALTY_INCREASE_FACTOR = 10.0;

          SquaredPenalty(ObjectiveFunction& f,
                         ObjectiveGradient& fGradient,
                         ConstraintFunction& g,
                         ConstraintGradient& gGradient,
                         ConstraintFunction& h,
                         ConstraintGradient& hGradient,
                         size_t maxItCount = DEFAULT_N,
                         float_t xTolerance = DEFAULT_X_TOLERANCE,
                         float_t constraintTolerance =
                           DEFAULT_CONSTRAINT_TOLERANCE,
                         float_t penaltyStartValue =
                           DEFAULT_PENALTY_START_VALUE,
                         float_t penaltyIncreaseFactor =
                           DEFAULT_PENALTY_INCREASE_FACTOR);

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

          /**
           * @return equality constraint function gradient
           */
          ConstraintGradient& getEqualityConstraintGradient() const;

          float_t getXTolerance() const;
          void setXTolerance(float_t xTolerance);

          float_t getConstraintTolerance() const;
          void setConstraintTolerance(float_t constraintTolerance);

          float_t getPenaltyStartValue() const;
          void setPenaltyStartValue(float_t penaltyStartValue);

          float_t getPenaltyIncreaseFactor() const;
          void setPenaltyIncreaseFactor(float_t penaltyIncreaseFactor);

        protected:
          /// objective function gradient
          ObjectiveGradient& fGradient;
          /// inequality constraint function gradient
          ConstraintGradient& gGradient;
          /// equality constraint function gradient
          ConstraintGradient& hGradient;
          float_t theta;
          float_t epsilon;
          float_t mu0;
          float_t rhoMuPlus;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_SQUAREDPENALTY_HPP */

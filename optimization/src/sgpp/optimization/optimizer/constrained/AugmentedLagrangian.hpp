// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_AUGMENTEDLAGRANGIAN_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_AUGMENTEDLAGRANGIAN_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/constrained/ConstrainedOptimizer.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionGradient.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Augmented Lagrangian method for constrained optimization.
       */
      class AugmentedLagrangian : public ConstrainedOptimizer {
        public:
          /// default point tolerance
          static constexpr float_t DEFAULT_X_TOLERANCE = 1e-6;
          /// default constraint tolerance
          static constexpr float_t DEFAULT_CONSTRAINT_TOLERANCE = 1e-6;
          /// default penalty start value
          static constexpr float_t DEFAULT_PENALTY_START_VALUE = 1.0;
          /// default penalty increase factor
          static constexpr float_t DEFAULT_PENALTY_INCREASE_FACTOR = 1.0;

          /**
           * Constructor.
           *
           * @param f                     objective function
           * @param fGradient             objective function gradient
           * @param g                     inequality constraint
           * @param gGradient             inequality constraint gradient
           * @param h                     equality constraint
           * @param hGradient             equality constraint gradient
           * @param maxItCount            maximal number of
           *                              function evaluations
           * @param xTolerance            point tolerance
           * @param constraintTolerance   constraint tolerance
           * @param penaltyStartValue     penalty start value
           * @param penaltyIncreaseFactor penalty increase factor
           */
          AugmentedLagrangian(ScalarFunction& f,
                              ScalarFunctionGradient& fGradient,
                              VectorFunction& g,
                              VectorFunctionGradient& gGradient,
                              VectorFunction& h,
                              VectorFunctionGradient& hGradient,
                              size_t maxItCount = DEFAULT_N,
                              float_t xTolerance = DEFAULT_X_TOLERANCE,
                              float_t constraintTolerance =
                                DEFAULT_CONSTRAINT_TOLERANCE,
                              float_t penaltyStartValue =
                                DEFAULT_PENALTY_START_VALUE,
                              float_t penaltyIncreaseFactor =
                                DEFAULT_PENALTY_INCREASE_FACTOR);

          void optimize();

          /**
           * Try to find a feasible initial point by solving an auxiliary
           * problem.
           * However, it is not guaranteed that the method works, i.e.,
           * it should be checked afterwards that the returned point is
           * actually feasible.
           *
           * @return feasible point in \f$[0, 1]^d\f$
           */
          base::DataVector findFeasiblePoint() const;

          /**
           * @return objective function gradient
           */
          ScalarFunctionGradient& getObjectiveGradient() const;

          /**
           * @return inequality constraint function gradient
           */
          VectorFunctionGradient& getInequalityConstraintGradient() const;

          /**
           * @return equality constraint function gradient
           */
          VectorFunctionGradient& getEqualityConstraintGradient() const;

          /**
           * @return point tolerance
           */
          float_t getXTolerance() const;

          /**
           * @param xTolerance point tolerance
           */
          void setXTolerance(float_t xTolerance);

          /**
           * @return constraint tolerance
           */
          float_t getConstraintTolerance() const;

          /**
           * @param constraintTolerance constraint tolerance
           */
          void setConstraintTolerance(float_t constraintTolerance);

          /**
           * @return penalty start value
           */
          float_t getPenaltyStartValue() const;

          /**
           * @param penaltyStartValue penalty start value
           */
          void setPenaltyStartValue(float_t penaltyStartValue);

          /**
           * @return penalty increase factor
           */
          float_t getPenaltyIncreaseFactor() const;

          /**
           * @param penaltyIncreaseFactor penalty increase factor
           */
          void setPenaltyIncreaseFactor(float_t penaltyIncreaseFactor);

        protected:
          /// objective function gradient
          ScalarFunctionGradient& fGradient;
          /// inequality constraint function gradient
          VectorFunctionGradient& gGradient;
          /// equality constraint function gradient
          VectorFunctionGradient& hGradient;
          /// point tolerance
          float_t theta;
          /// constraint tolerance
          float_t epsilon;
          /// penalty start value
          float_t mu0;
          /// penalty increase factor
          float_t rhoMuPlus;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_AUGMENTEDLAGRANGIAN_HPP */

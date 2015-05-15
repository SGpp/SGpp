// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_ADAPTIVENEWTON_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_ADAPTIVENEWTON_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/ObjectiveHessian.hpp>
#include <sgpp/optimization/sle/solver/GaussianElimination.hpp>
#include <sgpp/optimization/sle/solver/SLESolver.hpp>

#include <cstddef>
#include <memory>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      class AdaptiveNewton : public UnconstrainedOptimizer {
        public:
          static constexpr float_t DEFAULT_TOLERANCE = 1e-6;
          static constexpr float_t DEFAULT_STEP_SIZE_INCREASE_FACTOR = 1.2;
          static constexpr float_t DEFAULT_STEP_SIZE_DECREASE_FACTOR = 0.5;
          static constexpr float_t DEFAULT_DAMPING_INCREASE_FACTOR = 1.0;
          static constexpr float_t DEFAULT_DAMPING_DECREASE_FACTOR = 0.5;
          static constexpr float_t DEFAULT_LINE_SEARCH_ACCURACY = 0.01;

          AdaptiveNewton(ObjectiveFunction& f,
                         ObjectiveHessian& fHessian,
                         size_t maxItCount = DEFAULT_N,
                         float_t tolerance = DEFAULT_TOLERANCE,
                         float_t stepSizeIncreaseFactor = DEFAULT_STEP_SIZE_INCREASE_FACTOR,
                         float_t stepSizeDecreaseFactor = DEFAULT_STEP_SIZE_DECREASE_FACTOR,
                         float_t dampingIncreaseFactor = DEFAULT_DAMPING_INCREASE_FACTOR,
                         float_t dampingDecreaseFactor = DEFAULT_DAMPING_DECREASE_FACTOR,
                         float_t lineSearchAccuracy = DEFAULT_LINE_SEARCH_ACCURACY);

          AdaptiveNewton(ObjectiveFunction& f,
                         ObjectiveHessian& fHessian,
                         size_t maxItCount,
                         float_t tolerance,
                         float_t stepSizeIncreaseFactor,
                         float_t stepSizeDecreaseFactor,
                         float_t dampingIncreaseFactor,
                         float_t dampingDecreaseFactor,
                         float_t lineSearchAccuracy,
                         const sle_solver::SLESolver& sleSolver);

          /**
           * @param[out] xOpt optimal point
           * @return          optimal objective function value
           */
          float_t optimize(base::DataVector& xOpt);

          /**
           * @return objective function Hessian
           */
          ObjectiveHessian& getObjectiveHessian() const;

          float_t getTolerance() const;
          void setTolerance(float_t tolerance);

          float_t getStepSizeIncreaseFactor() const;
          void setStepSizeIncreaseFactor(float_t stepSizeIncreaseFactor);

          float_t getStepSizeDecreaseFactor() const;
          void setStepSizeDecreaseFactor(float_t stepSizeDecreaseFactor);

          float_t getDampingIncreaseFactor() const;
          void setDampingIncreaseFactor(float_t dampingIncreaseFactor);

          float_t getDampingDecreaseFactor() const;
          void setDampingDecreaseFactor(float_t dampingDecreaseFactor);

          float_t getLineSearchAccuracy() const;
          void setLineSearchAccuracy(float_t lineSearchAccuracy);

        protected:
          /// objective function Hessian
          ObjectiveHessian& fHessian;
          float_t theta;
          float_t rhoAlphaPlus;
          float_t rhoAlphaMinus;
          float_t rhoLambdaPlus;
          float_t rhoLambdaMinus;
          float_t rhoLs;
          /// default linear solver
          const sle_solver::GaussianElimination defaultSleSolver;
          /// linear solver
          const sle_solver::SLESolver& sleSolver;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_ADAPTIVENEWTON_HPP */

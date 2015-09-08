// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_ADAPTIVENEWTON_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_ADAPTIVENEWTON_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/optimization/sle/solver/GaussianElimination.hpp>
#include <sgpp/optimization/sle/solver/SLESolver.hpp>

#include <cstddef>
#include <memory>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Newton method with adaptive step size.
       */
      class AdaptiveNewton : public UnconstrainedOptimizer {
        public:
          /// default tolerance
          static constexpr float_t DEFAULT_TOLERANCE = 1e-6;
          /// default step size increase factor
          static constexpr float_t DEFAULT_STEP_SIZE_INCREASE_FACTOR = 1.2;
          /// default step size decrease factor
          static constexpr float_t DEFAULT_STEP_SIZE_DECREASE_FACTOR = 0.5;
          /// default damping increase factor
          static constexpr float_t DEFAULT_DAMPING_INCREASE_FACTOR = 1.0;
          /// default damping decrease factor
          static constexpr float_t DEFAULT_DAMPING_DECREASE_FACTOR = 0.5;
          /// default line search accuracy
          static constexpr float_t DEFAULT_LINE_SEARCH_ACCURACY = 0.01;

          /**
           * Constructor.
           * By default, GaussianElimination is used to solve the
           * linear systems.
           *
           * @param f                         objective function
           * @param fHessian                  objective function Hessian
           * @param maxItCount                maximal number of
           *                                  function evaluations
           * @param tolerance                 tolerance
           * @param stepSizeIncreaseFactor    step size increase factor
           * @param stepSizeDecreaseFactor    step size decrease factor
           * @param dampingIncreaseFactor     damping increase factor
           * @param dampingDecreaseFactor     damping decrease factor
           * @param lineSearchAccuracy        line search accuracy
           */
          AdaptiveNewton(ScalarFunction& f,
                         ScalarFunctionHessian& fHessian,
                         size_t maxItCount = DEFAULT_N,
                         float_t tolerance = DEFAULT_TOLERANCE,
                         float_t stepSizeIncreaseFactor = DEFAULT_STEP_SIZE_INCREASE_FACTOR,
                         float_t stepSizeDecreaseFactor = DEFAULT_STEP_SIZE_DECREASE_FACTOR,
                         float_t dampingIncreaseFactor = DEFAULT_DAMPING_INCREASE_FACTOR,
                         float_t dampingDecreaseFactor = DEFAULT_DAMPING_DECREASE_FACTOR,
                         float_t lineSearchAccuracy = DEFAULT_LINE_SEARCH_ACCURACY);

          /**
           * Constructor.
           * Do not destruct the solver before this object!
           *
           * @param f                         objective function
           * @param fHessian                  objective function Hessian
           * @param maxItCount                maximal number of
           *                                  function evaluations
           * @param tolerance                 tolerance
           * @param stepSizeIncreaseFactor    step size increase factor
           * @param stepSizeDecreaseFactor    step size decrease factor
           * @param dampingIncreaseFactor     damping increase factor
           * @param dampingDecreaseFactor     damping decrease factor
           * @param lineSearchAccuracy        line search accuracy
           * @param sleSolver                 reference to linear solver for
           *                                  solving the linear systems
           *                                  (Hessian as coefficient matrix)
           */
          AdaptiveNewton(ScalarFunction& f,
                         ScalarFunctionHessian& fHessian,
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
          ScalarFunctionHessian& getObjectiveHessian() const;

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
           * @return damping increase factor
           */
          float_t getDampingIncreaseFactor() const;

          /**
           * @param dampingIncreaseFactor damping increase factor
           */
          void setDampingIncreaseFactor(float_t dampingIncreaseFactor);

          /**
           * @return damping decrease factor
           */
          float_t getDampingDecreaseFactor() const;

          /**
           * @param dampingDecreaseFactor damping decrease factor
           */
          void setDampingDecreaseFactor(float_t dampingDecreaseFactor);

          /**
           * @return line search accuracy
           */
          float_t getLineSearchAccuracy() const;

          /**
           * @param lineSearchAccuracy line search accuracy
           */
          void setLineSearchAccuracy(float_t lineSearchAccuracy);

        protected:
          /// objective function Hessian
          ScalarFunctionHessian& fHessian;
          /// tolerance
          float_t theta;
          /// step size increase factor
          float_t rhoAlphaPlus;
          /// step size decrease factor
          float_t rhoAlphaMinus;
          /// damping increase factor
          float_t rhoLambdaPlus;
          /// damping decrease factor
          float_t rhoLambdaMinus;
          /// line search accuracy
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

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/sle/solver/GaussianElimination.hpp>
#include <sgpp/base/tools/sle/solver/SLESolver.hpp>
#include <sgpp/base/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

#include <cstddef>
#include <memory>

namespace sgpp {
namespace optimization {
namespace optimizer {

/**
 * Newton method with adaptive step size.
 */
class AdaptiveNewton : public UnconstrainedOptimizer {
 public:
  /// default tolerance
  static constexpr double DEFAULT_TOLERANCE = 1e-6;
  /// default step size increase factor
  static constexpr double DEFAULT_STEP_SIZE_INCREASE_FACTOR = 1.2;
  /// default step size decrease factor
  static constexpr double DEFAULT_STEP_SIZE_DECREASE_FACTOR = 0.5;
  /// default damping increase factor
  static constexpr double DEFAULT_DAMPING_INCREASE_FACTOR = 1.0;
  /// default damping decrease factor
  static constexpr double DEFAULT_DAMPING_DECREASE_FACTOR = 0.5;
  /// default line search accuracy
  static constexpr double DEFAULT_LINE_SEARCH_ACCURACY = 0.01;

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
  AdaptiveNewton(const base::ScalarFunction& f, const base::ScalarFunctionHessian& fHessian,
                 size_t maxItCount = DEFAULT_N, double tolerance = DEFAULT_TOLERANCE,
                 double stepSizeIncreaseFactor = DEFAULT_STEP_SIZE_INCREASE_FACTOR,
                 double stepSizeDecreaseFactor = DEFAULT_STEP_SIZE_DECREASE_FACTOR,
                 double dampingIncreaseFactor = DEFAULT_DAMPING_INCREASE_FACTOR,
                 double dampingDecreaseFactor = DEFAULT_DAMPING_DECREASE_FACTOR,
                 double lineSearchAccuracy = DEFAULT_LINE_SEARCH_ACCURACY);

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
  AdaptiveNewton(const base::ScalarFunction& f, const base::ScalarFunctionHessian& fHessian,
                 size_t maxItCount, double tolerance, double stepSizeIncreaseFactor,
                 double stepSizeDecreaseFactor, double dampingIncreaseFactor,
                 double dampingDecreaseFactor, double lineSearchAccuracy,
                 const base::sle_solver::SLESolver& sleSolver);

  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  AdaptiveNewton(const AdaptiveNewton& other);

  /**
   * Destructor.
   */
  ~AdaptiveNewton() override;

  void optimize() override;

  /**
   * @return tolerance
   */
  double getTolerance() const;

  /**
   * @param tolerance tolerance
   */
  void setTolerance(double tolerance);

  /**
   * @return step size increase factor
   */
  double getStepSizeIncreaseFactor() const;

  /**
   * @param stepSizeIncreaseFactor step size increase factor
   */
  void setStepSizeIncreaseFactor(double stepSizeIncreaseFactor);

  /**
   * @return step size decrease factor
   */
  double getStepSizeDecreaseFactor() const;

  /**
   * @param stepSizeDecreaseFactor step size decrease factor
   */
  void setStepSizeDecreaseFactor(double stepSizeDecreaseFactor);

  /**
   * @return damping increase factor
   */
  double getDampingIncreaseFactor() const;

  /**
   * @param dampingIncreaseFactor damping increase factor
   */
  void setDampingIncreaseFactor(double dampingIncreaseFactor);

  /**
   * @return damping decrease factor
   */
  double getDampingDecreaseFactor() const;

  /**
   * @param dampingDecreaseFactor damping decrease factor
   */
  void setDampingDecreaseFactor(double dampingDecreaseFactor);

  /**
   * @return line search accuracy
   */
  double getLineSearchAccuracy() const;

  /**
   * @param lineSearchAccuracy line search accuracy
   */
  void setLineSearchAccuracy(double lineSearchAccuracy);

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const override;

 protected:
  /// tolerance
  double theta;
  /// step size increase factor
  double rhoAlphaPlus;
  /// step size decrease factor
  double rhoAlphaMinus;
  /// damping increase factor
  double rhoLambdaPlus;
  /// damping decrease factor
  double rhoLambdaMinus;
  /// line search accuracy
  double rhoLs;
  /// default linear solver
  const base::sle_solver::GaussianElimination defaultSleSolver;
  /// linear solver
  const base::sle_solver::SLESolver& sleSolver;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

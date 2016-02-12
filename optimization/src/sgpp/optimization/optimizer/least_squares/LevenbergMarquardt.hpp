// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_LEAST_SQUARES_LEVENBERGMARQUARDT_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_LEAST_SQUARES_LEVENBERGMARQUARDT_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/least_squares/LeastSquaresOptimizer.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionGradient.hpp>
#include <sgpp/optimization/sle/solver/GaussianElimination.hpp>

namespace SGPP {
namespace optimization {
namespace optimizer {

/**
 * Levenberg-Marquardt algorithm for least squares optimization.
 */
class LevenbergMarquardt : public LeastSquaresOptimizer {
 public:
  /// default tolerance
  static constexpr float_t DEFAULT_TOLERANCE = 1e-6;
  /// default initial damping
  static constexpr float_t DEFAULT_INITIAL_DAMPING = 1.0;
  /// default acceptance threshold
  static constexpr float_t DEFAULT_ACCEPTANCE_THRESHOLD = 0.3;
  /// default effectiveness threshold
  static constexpr float_t DEFAULT_EFFECTIVENESS_THRESHOLD = 0.9;

  /**
   * Constructor.
   * By default, GaussianElimination is used to solve the
   * linear systems.
   *
   * @param phi                     base function
   * @param phiGradient             Jacobian of phi
   * @param maxItCount              maximal number of
   *                                function evaluations
   * @param tolerance               tolerance
   * @param initialDamping          initial damping
   * @param acceptanceThreshold     acceptance threshold
   * @param effectivenessThreshold  effectiveness threshold
   */
  LevenbergMarquardt(VectorFunction& phi,
                     VectorFunctionGradient& phiGradient,
                     size_t maxItCount = DEFAULT_N,
                     float_t tolerance = DEFAULT_TOLERANCE,
                     float_t initialDamping = DEFAULT_INITIAL_DAMPING,
                     float_t acceptanceThreshold =
                       DEFAULT_ACCEPTANCE_THRESHOLD,
                     float_t effectivenessThreshold =
                       DEFAULT_EFFECTIVENESS_THRESHOLD);

  /**
   * Constructor.
   * Do not destruct the solver before this object!
   *
   * @param phi                     phi function
   * @param phiGradient             Jacobian of phi
   * @param maxItCount              maximal number of
   *                                function evaluations
   * @param tolerance               tolerance
   * @param initialDamping          initial damping
   * @param acceptanceThreshold     acceptance threshold
   * @param effectivenessThreshold  effectiveness threshold
   * @param sleSolver               reference to linear solver
   *                                for solving the linear systems
   */
  LevenbergMarquardt(VectorFunction& phi,
                     VectorFunctionGradient& phiGradient,
                     size_t maxItCount,
                     float_t tolerance,
                     float_t initialDamping,
                     float_t acceptanceThreshold,
                     float_t effectivenessThreshold,
                     const sle_solver::SLESolver& sleSolver);

  /**
   * Destructor.
   */
  ~LevenbergMarquardt() override;

  void optimize() override;

  /**
   * @return phi gradient
   */
  VectorFunctionGradient& getPhiGradient() const;

  /**
   * @return tolerance
   */
  float_t getTolerance() const;

  /**
   * @param tolerance tolerance
   */
  void setTolerance(float_t tolerance);

  /**
   * @return initial damping
   */
  float_t getInitialDamping() const;

  /**
   * @param initialDamping initial damping
   */
  void setInitialDamping(float_t initialDamping);

  /**
   * @return acceptanceThreshold
   */
  float_t getAcceptanceThreshold() const;

  /**
   * @param acceptanceThreshold acceptance threshold
   */
  void setAcceptanceThreshold(float_t acceptanceThreshold);

  /**
   * @return effectiveness threshold
   */
  float_t getEffectivenessThreshold() const;

  /**
   * @param effectivenessThreshold effectiveness threshold
   */
  void setEffectivenessThreshold(float_t effectivenessThreshold);

  /**
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(std::unique_ptr<LeastSquaresOptimizer>& clone) const
  override;

 protected:
  /// phi gradient
  VectorFunctionGradient& phiGradient;
  /// tolerance
  float_t tol;
  /// initial damping
  float_t mu0;
  /// acceptance threshold
  float_t beta0;
  /// effectiveness threshold
  float_t beta1;
  /// default linear solver
  const sle_solver::GaussianElimination defaultSleSolver;
  /// linear solver
  const sle_solver::SLESolver& sleSolver;
};

}
}
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_LEAST_SQUARES_LEVENBERGMARQUARDT_HPP */

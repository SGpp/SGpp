// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_LEAST_SQUARES_LEVENBERGMARQUARDT_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_LEAST_SQUARES_LEVENBERGMARQUARDT_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/sle/solver/GaussianElimination.hpp>
#include <sgpp/optimization/optimizer/least_squares/LeastSquaresOptimizer.hpp>
#include "../../../../../../base/src/sgpp/base/function/vector/VectorFunctionGradient.hpp"

namespace sgpp {
namespace optimization {
namespace optimizer {

/**
 * Levenberg-Marquardt algorithm for least squares optimization.
 */
class LevenbergMarquardt : public LeastSquaresOptimizer {
 public:
  /// default tolerance
  static constexpr double DEFAULT_TOLERANCE = 1e-6;
  /// default initial damping
  static constexpr double DEFAULT_INITIAL_DAMPING = 1.0;
  /// default acceptance threshold
  static constexpr double DEFAULT_ACCEPTANCE_THRESHOLD = 0.3;
  /// default effectiveness threshold
  static constexpr double DEFAULT_EFFECTIVENESS_THRESHOLD = 0.9;

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
  LevenbergMarquardt(const base::VectorFunction& phi,
                     const base::VectorFunctionGradient& phiGradient, size_t maxItCount = DEFAULT_N,
                     double tolerance = DEFAULT_TOLERANCE,
                     double initialDamping = DEFAULT_INITIAL_DAMPING,
                     double acceptanceThreshold = DEFAULT_ACCEPTANCE_THRESHOLD,
                     double effectivenessThreshold = DEFAULT_EFFECTIVENESS_THRESHOLD);

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
  LevenbergMarquardt(const base::VectorFunction& phi,
                     const base::VectorFunctionGradient& phiGradient, size_t maxItCount,
                     double tolerance, double initialDamping, double acceptanceThreshold,
                     double effectivenessThreshold, const base::sle_solver::SLESolver& sleSolver);

  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  LevenbergMarquardt(const LevenbergMarquardt& other);

  /**
   * Destructor.
   */
  ~LevenbergMarquardt() override;

  void optimize() override;

  /**
   * @return phi gradient
   */
  base::VectorFunctionGradient& getPhiGradient() const;

  /**
   * @return tolerance
   */
  double getTolerance() const;

  /**
   * @param tolerance tolerance
   */
  void setTolerance(double tolerance);

  /**
   * @return initial damping
   */
  double getInitialDamping() const;

  /**
   * @param initialDamping initial damping
   */
  void setInitialDamping(double initialDamping);

  /**
   * @return acceptanceThreshold
   */
  double getAcceptanceThreshold() const;

  /**
   * @param acceptanceThreshold acceptance threshold
   */
  void setAcceptanceThreshold(double acceptanceThreshold);

  /**
   * @return effectiveness threshold
   */
  double getEffectivenessThreshold() const;

  /**
   * @param effectivenessThreshold effectiveness threshold
   */
  void setEffectivenessThreshold(double effectivenessThreshold);

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<LeastSquaresOptimizer>& clone) const override;

 protected:
  /// phi gradient
  std::unique_ptr<base::VectorFunctionGradient> phiGradient;
  /// tolerance
  double tol;
  /// initial damping
  double mu0;
  /// acceptance threshold
  double beta0;
  /// effectiveness threshold
  double beta1;
  /// default linear solver
  const base::sle_solver::GaussianElimination defaultSleSolver;
  /// linear solver
  const base::sle_solver::SLESolver& sleSolver;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_LEAST_SQUARES_LEVENBERGMARQUARDT_HPP */

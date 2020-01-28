// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_BFGS_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_BFGS_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

namespace sgpp {
namespace optimization {
namespace optimizer {

/**
 * BFGS method for unconstrained optimization.
 */
class BFGS : public UnconstrainedOptimizer {
 public:
  /// default tolerance
  static constexpr double DEFAULT_TOLERANCE = 1e-6;
  /// default step size increase factor
  static constexpr double DEFAULT_STEP_SIZE_INCREASE_FACTOR = 1.2;
  /// default step size decrease factor
  static constexpr double DEFAULT_STEP_SIZE_DECREASE_FACTOR = 0.5;
  /// default line search accuracy
  static constexpr double DEFAULT_LINE_SEARCH_ACCURACY = 0.01;

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
  BFGS(const base::ScalarFunction& f, const base::ScalarFunctionGradient& fGradient,
       size_t maxItCount = DEFAULT_N, double tolerance = DEFAULT_TOLERANCE,
       double stepSizeIncreaseFactor = DEFAULT_STEP_SIZE_INCREASE_FACTOR,
       double stepSizeDecreaseFactor = DEFAULT_STEP_SIZE_DECREASE_FACTOR,
       double lineSearchAccuracy = DEFAULT_LINE_SEARCH_ACCURACY);

  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  BFGS(const BFGS& other);

  /**
   * Destructor.
   */
  ~BFGS() override;

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
  /// line search accuracy
  double rhoLs;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_BFGS_HPP */

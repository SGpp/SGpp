// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_RPROP_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_RPROP_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

namespace sgpp {
namespace optimization {
namespace optimizer {

/**
 * Rprop method for unconstrained optimization.
 */
class Rprop : public UnconstrainedOptimizer {
 public:
  /// default tolerance
  static constexpr double DEFAULT_TOLERANCE = 1e-6;
  /// default initial step size
  static constexpr double DEFAULT_INITIAL_STEP_SIZE = 0.01;
  /// default step size increase factor
  static constexpr double DEFAULT_STEP_SIZE_INCREASE_FACTOR = 1.2;
  /// default step size decrease factor
  static constexpr double DEFAULT_STEP_SIZE_DECREASE_FACTOR = 0.5;

  /**
   * Constructor.
   *
   * @param f                       objective function
   * @param fGradient               objective function gradient
   * @param maxItCount              maximal number of
   *                                function evaluations
   * @param tolerance               tolerance
   * @param initialStepSize         initial step size
   * @param stepSizeIncreaseFactor  step size increase factor
   * @param stepSizeDecreaseFactor  step size decrease factor
   */
  Rprop(const base::ScalarFunction& f, const base::ScalarFunctionGradient& fGradient,
        size_t maxItCount = DEFAULT_N, double tolerance = DEFAULT_TOLERANCE,
        double initialStepSize = DEFAULT_INITIAL_STEP_SIZE,
        double stepSizeIncreaseFactor = DEFAULT_STEP_SIZE_INCREASE_FACTOR,
        double stepSizeDecreaseFactor = DEFAULT_STEP_SIZE_DECREASE_FACTOR);

  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  Rprop(const Rprop& other);

  /**
   * Destructor.
   */
  ~Rprop() override;

  void optimize() override;

  /**
   * @return objective function gradient
   */
  base::ScalarFunctionGradient& getObjectiveGradient() const;

  /**
   * @return tolerance
   */
  double getTolerance() const;

  /**
   * @param tolerance tolerance
   */
  void setTolerance(double tolerance);

  /**
   * @return initial step size
   */
  double getInitialStepSize() const;

  /**
   * @param initialStepSize initial step size
   */
  void setInitialStepSize(double initialStepSize);

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
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const override;

 protected:
  /// objective function gradient
  std::unique_ptr<base::ScalarFunctionGradient> fGradient;
  /// tolerance
  double theta;
  /// initial step size
  double initialAlpha;
  /// step size increase factor
  double rhoAlphaPlus;
  /// step size decrease factor
  double rhoAlphaMinus;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_RPROP_HPP */

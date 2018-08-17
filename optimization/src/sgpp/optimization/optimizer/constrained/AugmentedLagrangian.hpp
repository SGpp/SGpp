// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_AUGMENTEDLAGRANGIAN_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_AUGMENTEDLAGRANGIAN_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/constrained/ConstrainedOptimizer.hpp>

#include <vector>

namespace sgpp {
namespace optimization {
namespace optimizer {

/**
 * Augmented Lagrangian method for constrained optimization.
 */
class AugmentedLagrangian : public ConstrainedOptimizer {
 public:
  /// default point tolerance
  static constexpr double DEFAULT_X_TOLERANCE = 1e-6;
  /// default constraint tolerance
  static constexpr double DEFAULT_CONSTRAINT_TOLERANCE = 1e-6;
  /// default penalty start value
  static constexpr double DEFAULT_PENALTY_START_VALUE = 1.0;
  /// default penalty increase factor
  static constexpr double DEFAULT_PENALTY_INCREASE_FACTOR = 1.0;

  /**
   * Constructor with Nelder-Mead as optimization algorithm
   * (gradient-free).
   *
   * @param f                     objective function
   * @param g                     inequality constraint
   * @param h                     equality constraint
   * @param maxItCount            maximal number of
   *                              function evaluations
   * @param xTolerance            point tolerance
   * @param constraintTolerance   constraint tolerance
   * @param penaltyStartValue     penalty start value
   * @param penaltyIncreaseFactor penalty increase factor
   */
  AugmentedLagrangian(const ScalarFunction& f,
                      const VectorFunction& g,
                      const VectorFunction& h,
                      size_t maxItCount = DEFAULT_N,
                      double xTolerance = DEFAULT_X_TOLERANCE,
                      double constraintTolerance = DEFAULT_CONSTRAINT_TOLERANCE,
                      double penaltyStartValue = DEFAULT_PENALTY_START_VALUE,
                      double penaltyIncreaseFactor = DEFAULT_PENALTY_INCREASE_FACTOR);

  /**
   * Constructor with adaptive gradient descent as optimization algorithm
   * (gradient-based).
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
  AugmentedLagrangian(const ScalarFunction& f,
                      const ScalarFunctionGradient& fGradient,
                      const VectorFunction& g,
                      const VectorFunctionGradient& gGradient,
                      const VectorFunction& h,
                      const VectorFunctionGradient& hGradient,
                      size_t maxItCount = DEFAULT_N,
                      double xTolerance = DEFAULT_X_TOLERANCE,
                      double constraintTolerance = DEFAULT_CONSTRAINT_TOLERANCE,
                      double penaltyStartValue = DEFAULT_PENALTY_START_VALUE,
                      double penaltyIncreaseFactor = DEFAULT_PENALTY_INCREASE_FACTOR);

  /**
   * Constructor with custom unconstrained optimization algorithm
   * (gradient-free or gradient-based).
   *
   * @param unconstrainedOptimizer  unconstrained optimization algorithm
   * @param g                       inequality constraint
   * @param gGradient               inequality constraint gradient
   *                                (nullptr to omit)
   * @param h                       equality constraint
   * @param hGradient               equality constraint gradient
   *                                (nullptr to omit)
   * @param maxItCount              maximal number of
   *                                function evaluations
   * @param xTolerance              point tolerance
   * @param constraintTolerance     constraint tolerance
   * @param penaltyStartValue       penalty start value
   * @param penaltyIncreaseFactor   penalty increase factor
   */
  AugmentedLagrangian(const UnconstrainedOptimizer& unconstrainedOptimizer,
                      const VectorFunction& g,
                      const VectorFunctionGradient* gGradient,
                      const VectorFunction& h,
                      const VectorFunctionGradient* hGradient,
                      size_t maxItCount = DEFAULT_N,
                      double xTolerance = DEFAULT_X_TOLERANCE,
                      double constraintTolerance = DEFAULT_CONSTRAINT_TOLERANCE,
                      double penaltyStartValue = DEFAULT_PENALTY_START_VALUE,
                      double penaltyIncreaseFactor = DEFAULT_PENALTY_INCREASE_FACTOR);

  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  AugmentedLagrangian(const AugmentedLagrangian& other);

  /**
   * Destructor.
   */
  ~AugmentedLagrangian() override;

  void optimize() override;

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
   * @return point tolerance
   */
  double getXTolerance() const;

  /**
   * @param xTolerance point tolerance
   */
  void setXTolerance(double xTolerance);

  /**
   * @return constraint tolerance
   */
  double getConstraintTolerance() const;

  /**
   * @param constraintTolerance constraint tolerance
   */
  void setConstraintTolerance(double constraintTolerance);

  /**
   * @return penalty start value
   */
  double getPenaltyStartValue() const;

  /**
   * @param penaltyStartValue penalty start value
   */
  void setPenaltyStartValue(double penaltyStartValue);

  /**
   * @return penalty increase factor
   */
  double getPenaltyIncreaseFactor() const;

  /**
   * @param penaltyIncreaseFactor penalty increase factor
   */
  void setPenaltyIncreaseFactor(double penaltyIncreaseFactor);

  /**
   * @return tall matrix (d columns) in which the history of
   *         optimal points of the iterations are concatenated
   */
  const base::DataMatrix& getHistoryOfInnerIterationPoints() const;

  /**
   * @return vector in which the k-th entry indicates the number of
   *         inner iterations in the k-th (outer) iteration,
   *         empty vector on error
   */
  const std::vector<size_t>& getHistoryOfInnerIterationNumbers() const;

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const override;

 protected:
  /// point tolerance
  double theta;
  /// constraint tolerance
  double epsilon;
  /// penalty start value
  double mu0;
  /// penalty increase factor
  double rhoMuPlus;
  /// search history (inner iterations)
  base::DataMatrix xHistInner;
  /// search history (number of inner iterations)
  std::vector<size_t> kHistInner;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_AUGMENTEDLAGRANGIAN_HPP */

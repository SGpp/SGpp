// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_LOGBARRIER_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_LOGBARRIER_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/constrained/ConstrainedOptimizer.hpp>
#include <vector>
#include "../../../../../../base/src/sgpp/base/function/scalar/ScalarFunctionGradient.hpp"
#include "../../../../../../base/src/sgpp/base/function/vector/VectorFunctionGradient.hpp"

namespace sgpp {
namespace optimization {
namespace optimizer {

/**
 * Log Barrier method for constrained optimization.
 */
class LogBarrier : public ConstrainedOptimizer {
 public:
  /// default tolerance
  static constexpr double DEFAULT_TOLERANCE = 1e-6;
  /// default barrier start value
  static constexpr double DEFAULT_BARRIER_START_VALUE = 1.0;
  /// default barrier decrease factor
  static constexpr double DEFAULT_BARRIER_DECREASE_FACTOR = 0.5;

  /**
   * Constructor.
   *
   * @param f                     objective function
   * @param fGradient             objective function gradient
   * @param g                     inequality constraint
   * @param gGradient             inequality constraint gradient
   * @param maxItCount            maximal number of
   *                              function evaluations
   * @param tolerance             tolerance
   * @param barrierStartValue     barrier start value
   * @param barrierDecreaseFactor barrier decrease factor
   */
  LogBarrier(const base::ScalarFunction& f, const base::ScalarFunctionGradient& fGradient,
             const base::VectorFunction& g, const base::VectorFunctionGradient& gGradient,
             size_t maxItCount = DEFAULT_N, double tolerance = DEFAULT_TOLERANCE,
             double barrierStartValue = DEFAULT_BARRIER_START_VALUE,
             double barrierDecreaseFactor = DEFAULT_BARRIER_DECREASE_FACTOR);
  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  LogBarrier(const LogBarrier& other);

  /**
   * Destructor.
   */
  ~LogBarrier() override;

  void optimize() override;

  /**
   * @return objective function gradient
   */
  base::ScalarFunctionGradient& getObjectiveGradient() const;

  /**
   * @return inequality constraint function gradient
   */
  base::VectorFunctionGradient& getInequalityConstraintGradient() const;

  /**
   * @return tolerance
   */
  double getTolerance() const;

  /**
   * @param tolerance tolerance
   */
  void setTolerance(double tolerance);

  /**
   * @return barrier start value
   */
  double getBarrierStartValue() const;

  /**
   * @param barrierStartValue barrier start value
   */
  void setBarrierStartValue(double barrierStartValue);

  /**
   * @return barrier decrease factor
   */
  double getBarrierDecreaseFactor() const;

  /**
   * @param barrierDecreaseFactor barrier decrease factor
   */
  void setBarrierDecreaseFactor(double barrierDecreaseFactor);

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
  /// objective function gradient
  std::unique_ptr<base::ScalarFunctionGradient> fGradient;
  /// inequality constraint function gradient
  std::unique_ptr<base::VectorFunctionGradient> gGradient;
  /// tolerance
  double theta;
  /// barrier start value
  double mu0;
  /// barrier decrease factor
  double rhoMuMinus;
  /// search history (inner iterations)
  base::DataMatrix xHistInner;
  /// search history (number of inner iterations)
  std::vector<size_t> kHistInner;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_LOGBARRIER_HPP */

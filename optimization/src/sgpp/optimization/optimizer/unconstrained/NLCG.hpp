// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NLCG_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NLCG_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

namespace sgpp {
namespace optimization {
namespace optimizer {

/**
 * Gradient-based nonlinear conjugate gradient method.
 *
 * The method is restarted with the steepest descent direction
 * if the normalized absolute value of the inner product of
 * two successive gradients
 * exceeds a "restart threshold" \f$\alpha\f$.
 */
class NLCG : public UnconstrainedOptimizer {
 public:
  /// default beta (parameter for Armijo's rule)
  static constexpr double DEFAULT_BETA = 0.5;
  /// default gamma (parameter for Armijo's rule)
  static constexpr double DEFAULT_GAMMA = 1e-2;
  /// default tolerance (parameter for Armijo's rule)
  static constexpr double DEFAULT_TOLERANCE = 1e-8;
  /// default epsilon (parameter for Armijo's rule)
  static constexpr double DEFAULT_EPSILON = 1e-18;
  /// default restart threshold
  static constexpr double DEFAULT_RESTART_THRESHOLD = 0.1;

  /**
   * Constructor.
   *
   * @param f                 objective function
   * @param fGradient         objective function gradient
   * @param maxItCount        maximal number of iterations
   * @param beta              beta (parameter for Armijo's rule)
   * @param gamma             gamma (parameter for Armijo's rule)
   * @param tolerance         tolerance (parameter for Armijo's rule)
   * @param epsilon           epsilon (parameter for Armijo's rule)
   * @param restartThreshold  restart threshold
   */
  NLCG(const base::ScalarFunction& f, const base::ScalarFunctionGradient& fGradient,
       size_t maxItCount = DEFAULT_N, double beta = DEFAULT_BETA, double gamma = DEFAULT_GAMMA,
       double tolerance = DEFAULT_TOLERANCE, double epsilon = DEFAULT_EPSILON,
       double restartThreshold = DEFAULT_RESTART_THRESHOLD);

  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  NLCG(const NLCG& other);

  /**
   * Destructor.
   */
  ~NLCG() override;

  void optimize() override;

  /**
   * @return objective function gradient
   */
  base::ScalarFunctionGradient& getObjectiveGradient() const;

  /**
   * @return              beta (parameter for Armijo's rule)
   */
  double getBeta() const;

  /**
   * @param beta          beta (parameter for Armijo's rule)
   */
  void setBeta(double beta);

  /**
   * @return              gamma (parameter for Armijo's rule)
   */
  double getGamma() const;

  /**
   * @param gamma         gamma (parameter for Armijo's rule)
   */
  void setGamma(double gamma);

  /**
   * @return              tolerance (parameter for Armijo's rule)
   */
  double getTolerance() const;

  /**
   * @param tolerance     tolerance (parameter for Armijo's rule)
   */
  void setTolerance(double tolerance);

  /**
   * @return              epsilon (parameter for Armijo's rule)
   */
  double getEpsilon() const;

  /**
   * @param epsilon       epsilon (parameter for Armijo's rule)
   */
  void setEpsilon(double epsilon);

  /**
   * @return                  restart threshold
   */
  double getRestartThreshold() const;

  /**
   * @param restartThreshold  restart threshold
   */
  void setRestartThreshold(double restartThreshold);

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const override;

 protected:
  /// objective function gradient
  std::unique_ptr<base::ScalarFunctionGradient> fGradient;
  /// beta (parameter for Armijo's rule)
  double beta;
  /// gamma (parameter for Armijo's rule)
  double gamma;
  /// tolerance (parameter for Armijo's rule)
  double tol;
  /// epsilon (parameter for Armijo's rule)
  double eps;
  /// restart threshold
  double alpha;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NLCG_HPP */

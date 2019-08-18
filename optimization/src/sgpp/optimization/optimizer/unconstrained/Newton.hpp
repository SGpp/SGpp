// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NEWTON_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NEWTON_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/base/tools/sle/solver/GaussianElimination.hpp>
#include <sgpp/base/tools/sle/solver/SLESolver.hpp>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

#include <cstddef>
#include <memory>

namespace sgpp {
namespace optimization {
namespace optimizer {

/**
 * Gradient-based nonlinear conjugate gradient method.
 *
 * The method is restarted with the steepest descent direction
 * if the inner product of negated gradient and search direction is
 * not big enough (criterion depending on three parameters).
 */
class Newton : public UnconstrainedOptimizer {
 public:
  /// default beta (parameter for Armijo's rule)
  static constexpr double DEFAULT_BETA = 0.5;
  /// default gamma (parameter for Armijo's rule)
  static constexpr double DEFAULT_GAMMA = 1e-2;
  /// default tolerance (parameter for Armijo's rule)
  static constexpr double DEFAULT_TOLERANCE = 1e-8;
  /// default epsilon (parameter for Armijo's rule)
  static constexpr double DEFAULT_EPSILON = 1e-18;
  /// default steepest descent restart parameter 1
  static constexpr double DEFAULT_ALPHA1 = 1e-6;
  /// default steepest descent restart parameter 2
  static constexpr double DEFAULT_ALPHA2 = 1e-6;
  /// default steepest descent restart exponent
  static constexpr double DEFAULT_P = 0.1;

  /**
   * Constructor.
   * By default, GaussianElimination is used to solve the
   * linear systems.
   *
   * @param f                 objective function
   * @param fHessian          objective function Hessian
   * @param maxItCount        maximal number of iterations
   * @param beta              beta (parameter for Armijo's rule)
   * @param gamma             gamma (parameter for Armijo's rule)
   * @param tolerance         tolerance (parameter for Armijo's rule)
   * @param epsilon           epsilon (parameter for Armijo's rule)
   * @param alpha1            steepest descent restart parameter 1
   * @param alpha2            steepest descent restart parameter 2
   * @param p                 steepest descent restart exponent
   */
  Newton(const base::ScalarFunction& f, const base::ScalarFunctionHessian& fHessian,
         size_t maxItCount = DEFAULT_N, double beta = DEFAULT_BETA, double gamma = DEFAULT_GAMMA,
         double tolerance = DEFAULT_TOLERANCE, double epsilon = DEFAULT_EPSILON,
         double alpha1 = DEFAULT_ALPHA1, double alpha2 = DEFAULT_ALPHA2, double p = DEFAULT_P);

  /**
   * Constructor.
   * Do not destruct the solver before this object!
   *
   * @param f                 objective function
   * @param fHessian          objective function Hessian
   * @param maxItCount        maximal number of iterations
   * @param beta              beta (parameter for Armijo's rule)
   * @param gamma             gamma (parameter for Armijo's rule)
   * @param tolerance         tolerance (parameter for Armijo's rule)
   * @param epsilon           epsilon (parameter for Armijo's rule)
   * @param alpha1            steepest descent restart parameter 1
   * @param alpha2            steepest descent restart parameter 2
   * @param p                 steepest descent restart exponent
   * @param sleSolver         reference to linear solver for solving
   *                          the linear systems
   *                          (Hessian as coefficient matrix)
   */
  Newton(const base::ScalarFunction& f, const base::ScalarFunctionHessian& fHessian,
         size_t maxItCount, double beta, double gamma, double tolerance, double epsilon,
         double alpha1, double alpha2, double p, const base::sle_solver::SLESolver& sleSolver);

  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  Newton(const Newton& other);

  /**
   * Destructor.
   */
  ~Newton() override;

  void optimize() override;

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
   * @return          steepest descent restart parameter 1
   */
  double getAlpha1() const;

  /**
   * @param alpha1    steepest descent restart parameter 1
   */
  void setAlpha1(double alpha1);

  /**
   * @return          steepest descent restart parameter 2
   */
  double getAlpha2() const;

  /**
   * @param alpha2    steepest descent restart parameter 2
   */
  void setAlpha2(double alpha2);

  /**
   * @return          steepest descent restart exponent
   */
  double getP() const;

  /**
   * @param p         steepest descent restart exponent
   */
  void setP(double p);

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const override;

 protected:
  /// beta (parameter for Armijo's rule)
  double beta;
  /// gamma (parameter for Armijo's rule)
  double gamma;
  /// tolerance (parameter for Armijo's rule)
  double tol;
  /// epsilon (parameter for Armijo's rule)
  double eps;
  /// steepest descent restart parameter 1
  double alpha1;
  /// steepest descent restart parameter 2
  double alpha2;
  /// steepest descent restart exponent
  double p;
  /// default linear solver
  const base::sle_solver::GaussianElimination defaultSleSolver;
  /// linear solver
  const base::sle_solver::SLESolver& sleSolver;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NEWTON_HPP */

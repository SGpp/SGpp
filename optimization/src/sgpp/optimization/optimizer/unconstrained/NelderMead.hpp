// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NELDERMEAD_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NELDERMEAD_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

namespace sgpp {
namespace optimization {
namespace optimizer {

/**
 * Gradient-free Nelder-Mead method.
 */
class NelderMead : public UnconstrainedOptimizer {
 public:
  /// default reflection coefficient
  static constexpr double DEFAULT_ALPHA = 1.0;
  /// default expansion coefficient
  static constexpr double DEFAULT_BETA = 2.0;
  /// default contraction coefficient
  static constexpr double DEFAULT_GAMMA = 0.5;
  /// default shrinking coefficient
  static constexpr double DEFAULT_DELTA = 0.5;
  /// default maximal number of function evaluations
  static const size_t DEFAULT_MAX_FCN_EVAL_COUNT = 1000;
  /// edge length of starting simplex
  static constexpr double STARTING_SIMPLEX_EDGE_LENGTH = 0.4;

  /**
   * Constructor.
   * The starting point is set to
   * \f$(0.5, \dotsc, 0.5)^{\mathrm{T}}\f$.
   *
   * @param f                     objective function
   * @param maxFcnEvalCount       maximal number of
   *                              function evaluations
   * @param alpha                 reflection coefficient
   * @param beta                  expansion coefficient
   * @param gamma                 contraction coefficient
   * @param delta                 shrinking coefficient
   */
  NelderMead(const base::ScalarFunction& f, size_t maxFcnEvalCount = DEFAULT_MAX_FCN_EVAL_COUNT,
             double alpha = DEFAULT_ALPHA, double beta = DEFAULT_BETA, double gamma = DEFAULT_GAMMA,
             double delta = DEFAULT_DELTA);

  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  NelderMead(const NelderMead& other);

  /**
   * Destructor.
   */
  ~NelderMead() override;

  void optimize() override;

  /**
   * @return          reflection coefficient
   */
  double getAlpha() const;

  /**
   * @param alpha     reflection coefficient
   */
  void setAlpha(double alpha);

  /**
   * @return          expansion coefficient
   */
  double getBeta() const;

  /**
   * @param beta      expansion coefficient
   */
  void setBeta(double beta);

  /**
   * @return          contraction coefficient
   */
  double getGamma() const;

  /**
   * @param gamma     contraction coefficient
   */
  void setGamma(double gamma);

  /**
   * @return          shrinking coefficient
   */
  double getDelta() const;

  /**
   * @param delta     shrinking coefficient
   */
  void setDelta(double delta);

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const override;

 protected:
  /// reflection coefficient
  double alpha;
  /// expansion coefficient
  double beta;
  /// contraction coefficient
  double gamma;
  /// shrinking coefficient
  double delta;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NELDERMEAD_HPP */

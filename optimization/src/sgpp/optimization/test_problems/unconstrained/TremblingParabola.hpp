// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_TREMBLINGPARABOLA_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_TREMBLINGPARABOLA_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Trembling parabola objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * \sum_{t=1}^d \left(\bar{x}_t^2/10 +
 * \frac{|\bar{x}_t|}{2}
 * \left(\frac{g^p(\bar{x}_t)}{g^p(0)} + 1\right)\right),\qquad
 * g^p(x) := \sum_{k \sim x} (-1)^k b^p(x - k + (p+1)/2)\f]
 * with \f$b^p\f$ being the B-spline with knots
 * \f$(0, 1, \dotsc, p + 1)\f$ (degree \f$p\f$) and
 * the sum running over all \f$k\f$ for which the summand doesn't vanish
 */
class TremblingParabolaObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   * @param p     smoothness of the function
   *              (\f$p\f$ times continuously differentiable)
   */
  TremblingParabolaObjective(size_t d, size_t p);

  /**
   * Destructor.
   */
  ~TremblingParabolaObjective() override;

  /**
   * @param x     point \f$\vec{x} \in [0, 1]^d\f$
   * @return      \f$f(\vec{x})\f$
   */
  double evalUndisplaced(const base::DataVector& x) override;

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<base::ScalarFunction>& clone) const override;

 protected:
  const size_t p;
  const base::SBsplineBase bSplineBasis;
  const double g0;

  inline double splineTrembling(double x) const;
};

/**
 * Trembling parabola unconstrained test problem.
 *
 * * Number of parameters: \f$d\f$
 * * Domain: \f$\bar{\vec{x}} \in [-4, 16]^d\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   \vec{0}\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   0\f$
 */
class TremblingParabola : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   * @param p     smoothness of the function
   *              (\f$p\f$ times continuously differentiable)
   */
  TremblingParabola(size_t d, size_t p);

  /**
   * Destructor.
   */
  ~TremblingParabola() override;

  /**
   * @return  objective function of the test problem
   */
  TestScalarFunction& getObjectiveFunction() override;

  /**
   * @param[out] x minimal point
   *               \f$\vec{x}_\opt \in [0, 1]^d\f$
   * @return       minimal function value
   *               \f$f(\vec{x}_\opt)\f$
   */
  double getOptimalPointUndisplaced(base::DataVector& x) override;

 protected:
  /// objective function
  TremblingParabolaObjective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_TREMBLINGPARABOLA_HPP */

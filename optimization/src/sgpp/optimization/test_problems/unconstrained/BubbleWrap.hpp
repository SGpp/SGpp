// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_BUBBLEWRAP_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_BUBBLEWRAP_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Bubble wrap objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * 1 - \prod_{t=1}^d (9/10 - |\bar{x}_t| + \cos(10\pi \bar{x}_t)/10)\f]
 */
class BubbleWrapObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  explicit BubbleWrapObjective(size_t d);

  /**
   * Destructor.
   */
  ~BubbleWrapObjective() override;

  /**
   * @param x     point \f$\vec{x} \in [0, 1]^d\f$
   * @return      \f$f(\vec{x})\f$
   */
  double evalUndisplaced(const base::DataVector& x) override;

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<base::ScalarFunction>& clone) const override;
};

/**
 * Bubble wrap unconstrained test problem.
 *
 * * Number of parameters: \f$d\f$
 * * Domain: \f$\bar{\vec{x}} \in d\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   \vec{0}\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   0\f$
 */
class BubbleWrap : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  explicit BubbleWrap(size_t d);

  /**
   * Destructor.
   */
  ~BubbleWrap() override;

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
  BubbleWrapObjective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_BUBBLEWRAP_HPP */

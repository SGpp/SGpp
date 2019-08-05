// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_ABSOLUTEVALUE_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_ABSOLUTEVALUE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Absolute value objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) := \sum_{t=1}^d |\bar{x}_t - 2^{-t}|\f]
 */
class AbsoluteValueObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  explicit AbsoluteValueObjective(size_t d);

  /**
   * Destructor.
   */
  ~AbsoluteValueObjective() override;

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
 * Absolute value unconstrained test problem.
 *
 * * Number of parameters: 2
 * * Domain: \f$\bar{\vec{x}} \in [0, 1]^d\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (2^{-t})_{t=1}^d\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   0\f$
 */
class AbsoluteValue : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  explicit AbsoluteValue(size_t d);

  /**
   * Destructor.
   */
  ~AbsoluteValue() override;

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
  AbsoluteValueObjective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_ABSOLUTEVALUE_HPP */

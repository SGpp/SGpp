// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_SHCB_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_SHCB_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * SHCB objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * \bar{x}_1^2 \left(4 - 2.1 \bar{x}_1^2 + \bar{x}_1^4/3\right) +
 * \bar{x}_1 \bar{x}_2 +
 * 4 \bar{x}_2^2 \left(\bar{x}_2^2 - 1\right)\f]
 */
class SHCBObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  SHCBObjective();

  /**
   * Destructor.
   */
  ~SHCBObjective() override;

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
 * SHCB unconstrained test problem.
 *
 * * Number of parameters: 2
 * * Domain: \f$\bar{\vec{x}} \in [-5, 5]^2\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} \in
 *   \{(0.08984247, -0.7126564), (-0.08984247, 0.7126564)\}\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -1.031628\f$
 */
class SHCB : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  SHCB();

  /**
   * Destructor.
   */
  ~SHCB() override;

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
  SHCBObjective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_SHCB_HPP */

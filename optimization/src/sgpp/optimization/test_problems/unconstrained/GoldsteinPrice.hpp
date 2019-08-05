// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_GOLDSTEINPRICE_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_GOLDSTEINPRICE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Goldstein-Price objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * 10^{-4} \cdot
 * (1 + (\bar{x}_1+\bar{x}_2+1)^2 (19 - 14\bar{x}_1 + 3\bar{x}_1^2 -
 * 14\bar{x}_2 + 6\bar{x}_1 \bar{x}_2 + 3\bar{x}_2^2))
 * \cdot (30 + (2\bar{x}_1 - 3\bar{x}_2)^2
 * (18 - 32\bar{x}_1 + 12\bar{x}_1^2 + 48\bar{x}_2 -
 * 36\bar{x}_1 \bar{x}_2 + 27\bar{x}_2^2))\f]
 */
class GoldsteinPriceObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  GoldsteinPriceObjective();

  /**
   * Destructor.
   */
  ~GoldsteinPriceObjective() override;

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
 * Goldstein-Price unconstrained test problem.
 *
 * * Number of parameters: 2
 * * Domain: \f$\bar{\vec{x}} \in [-2, 2]^2\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (0, -1)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   0.0003\f$
 */
class GoldsteinPrice : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  GoldsteinPrice();

  /**
   * Destructor.
   */
  ~GoldsteinPrice() override;

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
  GoldsteinPriceObjective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_GOLDSTEINPRICE_HPP */

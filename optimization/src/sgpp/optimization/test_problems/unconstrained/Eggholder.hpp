// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_EGGHOLDER_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_EGGHOLDER_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Eggholder objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x})} :=
 * -(\bar{x}_2 + 47) \sin\!\left(\sqrt{
 * \left|\bar{x}_1/2 + \bar{x}_2 + 47\right|}\right) -
 * \bar{x}_1 \sin\!\left(\sqrt{
 * \left|\bar{x}_1 - \bar{x}_2 - 47\right|}\right)\f]
 */
class EggholderObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  EggholderObjective();

  /**
   * Destructor.
   */
  ~EggholderObjective() override;

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
 * Eggholder unconstrained test problem.
 *
 * * Number of parameters: 2
 * * Domain: \f$\bar{\vec{x}} \in [-512, 512]^2\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (512, 404.2318)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -959.6407\f$
 *
 * The displacement is restricted because the minimal point lies on
 * the boundary of \f$[0, 1]^2\f$.
 */
class Eggholder : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  Eggholder();

  /**
   * Destructor.
   */
  ~Eggholder() override;

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
  /**
   * Sets the first displacement component to zero and
   * checks the resulting displacement.
   *
   * @return whether the current displacement is feasible
   */
  bool isDisplacementFeasible() override;

  /// objective function
  EggholderObjective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_EGGHOLDER_HPP */

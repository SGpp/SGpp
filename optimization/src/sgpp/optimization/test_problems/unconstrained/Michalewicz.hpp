// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_MICHALEWICZ_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_MICHALEWICZ_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Michalewicz objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * -\sin \bar{x}_1 \sin^{20}\!\left(\bar{x}_1^2/\pi\right) -
 * \sin \bar{x}_2 \sin^{20}\!\left(2\bar{x}_2^2/\pi\right)\f]
 */
class MichalewiczObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  MichalewiczObjective();

  /**
   * Destructor.
   */
  ~MichalewiczObjective() override;

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
 * Michalewicz unconstrained test problem.
 *
 * * Number of parameters: 2
 * * Domain: \f$\bar{\vec{x}} \in [0, 5]^2\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (2.202906, \pi/2)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -1.801303\f$
 */
class Michalewicz : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  Michalewicz();

  /**
   * Destructor.
   */
  ~Michalewicz() override;

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
  MichalewiczObjective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_MICHALEWICZ_HPP */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_MLADINEO_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_MLADINEO_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Mladineo objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * 1 + \frac{1}{2} \norm{\bar{\vec{x}}}_2^2 -
 * \cos(10 \ln(2\bar{x}_1)) \cos(10 \ln(3\bar{x}_2))\f]
 */
class MladineoObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  MladineoObjective();

  /**
   * Destructor.
   */
  ~MladineoObjective() override;

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
 * Mladineo unconstrained test problem.
 *
 * * Number of parameters: 2
 * * Domain: \f$\bar{\vec{x}} \in [0.01, 1]^2\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (0.01152704, 0.01440461)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   0.0001701830\f$
 */
class Mladineo : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  Mladineo();

  /**
   * Destructor.
   */
  ~Mladineo() override;

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
   * Checks if \f$\vec{d} \in [-0.01, 0] \times [-0.01, 0]\f$.
   *
   * @return whether the current displacement is feasible
   */
  bool isDisplacementFeasible() override;

  /// objective function
  MladineoObjective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_MLADINEO_HPP */

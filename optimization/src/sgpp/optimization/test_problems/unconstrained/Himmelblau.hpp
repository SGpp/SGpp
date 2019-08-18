// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_HIMMELBLAU_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_HIMMELBLAU_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Himmelblau objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * (\bar{x}_1^2 + \bar{x}_2 - 11)^2 + (\bar{x}_1 + \bar{x}_2^2 - 7)^2\f]
 */
class HimmelblauObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  HimmelblauObjective();

  /**
   * Destructor.
   */
  ~HimmelblauObjective() override;

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
 * Himmelblau unconstrained test problem.
 *
 * * Number of parameters: 2
 * * Domain: \f$\bar{\vec{x}} \in [-5, 5]^2\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} \in
 *   \{(3, 2), (-2.805118, 3.131313), (-3.779310, -3.283186),
 *   (3.584428, -1.848127)\}\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   0\f$
 */
class Himmelblau : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  Himmelblau();

  /**
   * Destructor.
   */
  ~Himmelblau() override;

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
  HimmelblauObjective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_HIMMELBLAU_HPP */

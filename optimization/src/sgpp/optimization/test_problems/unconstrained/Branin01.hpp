// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_BRANIN01_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_BRANIN01_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Branin01 objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * \left(\bar{x}_2 - 5.1 \bar{x}_1^2/(4\pi^2) +
 * 5 \bar{x}_1/\pi - 6\right)^2 +
 * 10 \left(1 - 1/(8\pi)\right) \cos \bar{x}_1 + 10\f]
 */
class Branin01Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  Branin01Objective();

  /**
   * Destructor.
   */
  ~Branin01Objective() override;

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
 * Branin01 unconstrained test problem.
 *
 * * Number of parameters: 2
 * * Domain: \f$\bar{\vec{x}} \in [-5, 10] \times [0, 15]\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} \in
 *   \{(-\pi, 491/40), (\pi, 91/40), (3\pi, 99/40)\}\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   5/(4\pi)\f$
 */
class Branin01 : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  Branin01();

  /**
   * Destructor.
   */
  ~Branin01() override;

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
  Branin01Objective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_BRANIN01_HPP */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_GRIEWANK_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_GRIEWANK_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Griewank objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * 1 + \frac{\norm{\bar{\vec{x}}}_2^2}{4000} -
 * \prod_{t=1}^d \cos\!\left(\frac{\bar{x}_t}{\sqrt{t}}\right)\f]
 */
class GriewankObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  explicit GriewankObjective(size_t d);

  /**
   * Destructor.
   */
  ~GriewankObjective() override;

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
 * Griewank unconstrained test problem.
 *
 * * Number of parameters: \f$d\f$
 * * Domain: \f$\bar{\vec{x}} \in [-600, 600]^d\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   \vec{0}\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   0\f$
 */
class Griewank : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  explicit Griewank(size_t d);

  /**
   * Destructor.
   */
  ~Griewank() override;

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
  GriewankObjective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_GRIEWANK_HPP */

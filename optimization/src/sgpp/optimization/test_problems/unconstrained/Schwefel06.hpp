// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_SCHWEFEL06_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_SCHWEFEL06_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Schwefel06 objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * \max(|\bar{x}_1 + 2 \bar{x}_2 - 7|, |2 \bar{x}_1 + \bar{x}_2 - 5|)\f]
 */
class Schwefel06Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  Schwefel06Objective();

  /**
   * Destructor.
   */
  ~Schwefel06Objective() override;

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
 * Schwefel06 unconstrained test problem.
 *
 * * Number of parameters: 2
 * * Domain: \f$\bar{\vec{x}} \in [-100, 100]^2\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} = (1, 3)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) = 0\f$
 */
class Schwefel06 : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  Schwefel06();

  /**
   * Destructor.
   */
  ~Schwefel06() override;

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
  Schwefel06Objective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_SCHWEFEL06_HPP */

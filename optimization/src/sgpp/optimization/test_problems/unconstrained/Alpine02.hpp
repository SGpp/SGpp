// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_ALPINE02_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_ALPINE02_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Alpine02 objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * -\prod_{t=1}^d \sqrt{\bar{x}_t} \sin(\bar{x}_t)\f]
 */
class Alpine02Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  explicit Alpine02Objective(size_t d);

  /**
   * Destructor.
   */
  ~Alpine02Objective() override;

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
 * Alpine02 unconstrained test problem.
 *
 * * Number of parameters: \f$d\f$
 * * Domain: \f$\bar{\vec{x}} \in [2, 10]^d\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   7.917052684666 \cdot \vec{1}\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -2.8081311800070^d\f$
 */
class Alpine02 : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  explicit Alpine02(size_t d);

  /**
   * Destructor.
   */
  ~Alpine02() override;

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
  Alpine02Objective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_ALPINE02_HPP */

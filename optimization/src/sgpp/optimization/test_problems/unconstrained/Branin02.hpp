// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_BRANIN02_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_BRANIN02_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

#include <cmath>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Branin02 objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * \left(\bar{x}_2 - 5.1 \bar{x}_1^2/(4\pi^2) +
 * 5 \bar{x}_1/\pi - 6\right)^2 +
 * 10 \left(1 - 1/(8\pi)\right) \cos \bar{x}_1 \cos \bar{x}_2 +
 * \ln(\bar{x}_1^2 + \bar{x}_2^2 + 1) + 10\f]
 */
class Branin02Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  Branin02Objective();

  /**
   * Destructor.
   */
  ~Branin02Objective() override;

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
 * Branin02 unconstrained test problem.
 *
 * * Number of parameters: 2
 * * Domain: \f$\bar{\vec{x}} \in [-5, 15]^2\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (-3.196988424804, 12.52625788532)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   5.559037320859\f$
 */
class Branin02 : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  Branin02();

  /**
   * Destructor.
   */
  ~Branin02() override;

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
  Branin02Objective f;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_BRANIN02_HPP */

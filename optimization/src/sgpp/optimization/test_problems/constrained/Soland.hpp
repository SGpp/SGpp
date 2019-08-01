// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_SOLAND_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_SOLAND_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Soland objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * -12 \bar{x}_1 - 7 \bar{x}_2 + \bar{x}_2^2\f]
 */
class SolandObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  SolandObjective();

  /**
   * Destructor.
   */
  ~SolandObjective() override;

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
 * Soland inequality constraint function.
 *
 * Definition: empty, i.e., no constraint
 */
class SolandInequalityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  SolandInequalityConstraint();

  /**
   * Destructor.
   */
  ~SolandInequalityConstraint() override;

  /**
   * @param       x       point \f$\vec{x} \in \mathbb{R}^d\f$
   * @param[out]  value   \f$\vec{f}(\vec{x})\f$
   */
  void evalUndisplaced(const base::DataVector& x, base::DataVector& value) override;

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<base::VectorFunction>& clone) const override;
};

/**
 * Soland equality constraint function.
 *
 * Definition:
 * \f[\bar{h}(\bar{\vec{x}}) :=
 * -2 \bar{x}_1^4 - \bar{x}_2 + 2\f]
 */
class SolandEqualityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  SolandEqualityConstraint();

  /**
   * Destructor.
   */
  ~SolandEqualityConstraint() override;

  /**
   * @param       x       point \f$\vec{x} \in \mathbb{R}^d\f$
   * @param[out]  value   \f$\vec{f}(\vec{x})\f$
   */
  void evalUndisplaced(const base::DataVector& x, base::DataVector& value) override;

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<base::VectorFunction>& clone) const override;
};

/**
 * Soland constrained test problem.
 *
 * * Number of parameters: 2
 * * Number of inequality constraints: 0
 * * Number of equality constraints: 1
 * * Domain: \f$\bar{\vec{x}} \in [0, 2] \times [0, 3]\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (0.7175362, 1.469842)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -16.73889\f$
 */
class Soland : public ConstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  Soland();

  /**
   * Destructor.
   */
  ~Soland() override;

  /**
   * @return  objective function of the test problem
   */
  TestScalarFunction& getObjectiveFunction() override;

  /**
   * @return  inequality function of the test problem
   */
  TestVectorFunction& getInequalityConstraintFunction() override;

  /**
   * @return  equality constraint function of the test problem
   */
  TestVectorFunction& getEqualityConstraintFunction() override;

  /**
   * @param[out] x minimal point
   *               \f$\vec{x}_\opt \in [0, 1]^d\f$
   * @return       minimal function value
   *               \f$f(\vec{x}_\opt)\f$
   */
  double getOptimalPointUndisplaced(base::DataVector& x) override;

 protected:
  /// objective function
  SolandObjective f;
  /// inequality constraint function
  SolandInequalityConstraint g;
  /// equality constraint function
  SolandEqualityConstraint h;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_SOLAND_HPP */

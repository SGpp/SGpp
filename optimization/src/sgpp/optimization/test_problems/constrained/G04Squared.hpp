// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G04SQUARED_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G04SQUARED_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * G04Squared objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * (5.3578547 \bar{x}_3^2 + 0.8356891 \bar{x}_1 \bar{x}_5 +
 * 37.293239 \bar{x}_1 - 10120)^2\f]
 */
class G04SquaredObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  G04SquaredObjective();

  /**
   * Destructor.
   */
  ~G04SquaredObjective() override;

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
 * G04Squared inequality constraint function.
 *
 * Definition:
 * \f[\bar{\vec{g}}(\bar{\vec{x}}) :=
 * \begin{pmatrix}
 *     85.334407 + 0.0056858 \bar{x}_2 \bar{x}_5 +
 *     0.0006262 \bar{x}_1 \bar{x}_4 - 0.0022053 \bar{x}_3 \bar{x}_5 -
 *     92\\
 *     -85.334407 - 0.0056858 \bar{x}_2 \bar{x}_5 -
 *     0.0006262 \bar{x}_1 \bar{x}_4 + 0.0022053 \bar{x}_3 \bar{x}_5\\
 *     80.51249 + 0.0071317 \bar{x}_2 \bar{x}_5 +
 *     0.0029955 \bar{x}_1 \bar{x}_2 + 0.0021813 \bar{x}_3^2 - 110\\
 *     -80.51249 - 0.0071317 \bar{x}_2 \bar{x}_5 -
 *     0.0029955 \bar{x}_1 \bar{x}_2 - 0.0021813 \bar{x}_3^2 + 90\\
 *     9.300961 + 0.0047026 \bar{x}_3 \bar{x}_5 +
 *     0.0012547 \bar{x}_1 \bar{x}_3 + 0.0019085 \bar{x}_3 \bar{x}_4 -
 *     25\\
 *     -9.300961 - 0.0047026 \bar{x}_3 \bar{x}_5 -
 *     0.0012547 \bar{x}_1 \bar{x}_3 - 0.0019085 \bar{x}_3 \bar{x}_4 +
 *     20
 * \end{pmatrix}\f]
 */
class G04SquaredInequalityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G04SquaredInequalityConstraint();

  /**
   * Destructor.
   */
  ~G04SquaredInequalityConstraint() override;

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
 * G04Squared equality constraint function.
 *
 * Definition: empty, i.e., no constraint
 */
class G04SquaredEqualityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G04SquaredEqualityConstraint();

  /**
   * Destructor.
   */
  ~G04SquaredEqualityConstraint() override;

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
 * G04Squared constrained test problem.
 *
 * * Number of parameters: 5
 * * Number of inequality constraints: 6
 * * Number of equality constraints: 0
 * * Domain: \f$\bar{\vec{x}} \in
 *   [78, 102] \times [33, 45] \times [27, 45]^3\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (78, 33, 29.99526, 45, 36.77581)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   43.590738\f$
 */
class G04Squared : public ConstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  G04Squared();

  /**
   * Destructor.
   */
  ~G04Squared() override;

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
  /**
   * Sets the first, second, and fourth displacement component to zero and
   * checks the resulting displacement.
   *
   * @return whether the current displacement is feasible
   */
  bool isDisplacementFeasible() override;

  /// objective function
  G04SquaredObjective f;
  /// inequality constraint function
  G04SquaredInequalityConstraint g;
  /// equality constraint function
  G04SquaredEqualityConstraint h;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G04SQUARED_HPP */

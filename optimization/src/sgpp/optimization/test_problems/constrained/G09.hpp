// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G09_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G09_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * G09 objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * (\bar{x}_1 - 10)^2 + 5 (\bar{x}_2 - 12)^2 + \bar{x}_3^4 +
 * 3 (\bar{x}_4 - 11)^2 + 10 \bar{x}_5^6 + 7 \bar{x}_6^2 + \bar{x}_7^4 -
 * 4 \bar{x}_6 \bar{x}_7 - 10 \bar{x}_6 - 8 \bar{x}_7\f]
 */
class G09Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  G09Objective();

  /**
   * Destructor.
   */
  ~G09Objective() override;

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
 * G09 inequality constraint function.
 *
 * Definition:
 * \f[\bar{\vec{g}}(\bar{\vec{x}}) :=
 * \begin{pmatrix}
 *     -127 + 2 \bar{x}_1^2 + 3 \bar{x}_2^4 + \bar{x}_3
 *     + 4 \bar{x}_4^2 + 5 \bar{x}_2\\
 *     -282 + 7 \bar{x}_1 + 3 \bar{x}_2 + 10 \bar{x}_3^2
 *     + \bar{x}_4 - \bar{x}_5\\
 *     -196 + 23 \bar{x}_1 + \bar{x}_2^2 + 6 \bar{x}_6^2 -
 *     8 \bar{x}_7\\
 *     4 \bar{x}_1^2 + \bar{x}_2^2 - 3 \bar{x}_1 \bar{x}_2 +
 *     2 \bar{x}_3^2 + 5 \bar{x}_6 - 11 \bar{x}_7
 * \end{pmatrix}\f]
 */
class G09InequalityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G09InequalityConstraint();

  /**
   * Destructor.
   */
  ~G09InequalityConstraint() override;

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
 * G09 equality constraint function.
 *
 * Definition: empty, i.e., no constraint
 */
class G09EqualityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G09EqualityConstraint();

  /**
   * Destructor.
   */
  ~G09EqualityConstraint() override;

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
 * G09 constrained test problem.
 *
 * * Number of parameters: 7
 * * Number of inequality constraints: 4
 * * Number of equality constraints: 0
 * * Domain: \f$\bar{\vec{x}} \in [-10, 10]^7\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (-2.330499, 1.951372, -0.4775414, 4.365726, 1.038131, 1.594227)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   680.6301\f$
 */
class G09 : public ConstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  G09();

  /**
   * Destructor.
   */
  ~G09() override;

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
  G09Objective f;
  /// inequality constraint function
  G09InequalityConstraint g;
  /// equality constraint function
  G09EqualityConstraint h;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G09_HPP */

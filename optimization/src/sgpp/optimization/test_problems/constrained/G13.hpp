// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G13_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G13_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * G13 objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * \exp\!\Big(\prod_{t=1}^5 \bar{x}_t\Big)\f]
 */
class G13Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  G13Objective();

  /**
   * Destructor.
   */
  ~G13Objective() override;

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
 * G13 inequality constraint function.
 *
 * Definition: empty, i.e., no constraint
 */
class G13InequalityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G13InequalityConstraint();

  /**
   * Destructor.
   */
  ~G13InequalityConstraint() override;

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
 * G13 equality constraint function.
 *
 * Definition:
 * \f[\bar{\vec{h}}(\bar{\vec{x}}) :=
 * \begin{pmatrix}
 *     -10 + \norm{\bar{\vec{x}}}_2^2\\
 *     \bar{x}_2 \bar{x}_3 - 5 \bar{x}_4 \bar{x}_5\\
 *     \bar{x}_1^3 + \bar{x}_2^3 + 1
 * \end{pmatrix}\f]
 */
class G13EqualityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G13EqualityConstraint();

  /**
   * Destructor.
   */
  ~G13EqualityConstraint() override;

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
 * G13 constrained test problem.
 *
 * * Number of parameters: 5
 * * Number of inequality constraints: 0
 * * Number of equality constraints: 3
 * * Domain: \f$\bar{\vec{x}} \in [-2.3, 2.3]^2 \times [-3.2, 3.2]^3\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (-1.717143, 1.595709, 1.827247, -0.7636413, -0.7636446)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   0.05394991\f$
 */
class G13 : public ConstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  G13();

  /**
   * Destructor.
   */
  ~G13() override;

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
  G13Objective f;
  /// inequality constraint function
  G13InequalityConstraint g;
  /// equality constraint function
  G13EqualityConstraint h;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G13_HPP */

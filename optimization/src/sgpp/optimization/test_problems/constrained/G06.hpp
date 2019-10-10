// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G06_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G06_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * G06 objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * (\bar{x}_1 - 10)^3 + (\bar{x}_2 - 20)^3\f]
 */
class G06Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  G06Objective();

  /**
   * Destructor.
   */
  ~G06Objective() override;

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
 * G06 inequality constraint function.
 *
 * Definition:
 * \f[\bar{\vec{g}}(\bar{\vec{x}}) :=
 * \begin{pmatrix}
 *     -(\bar{x}_1 - 5)^2 - (\bar{x}_2 - 5)^2 + 100\\
 *     (\bar{x}_1 - 6)^2 + (\bar{x}_2 - 5)^2 - 82.81
 * \end{pmatrix}\f]
 */
class G06InequalityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G06InequalityConstraint();

  /**
   * Destructor.
   */
  ~G06InequalityConstraint() override;

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
 * G06 equality constraint function.
 *
 * Definition: empty, i.e., no constraint
 */
class G06EqualityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G06EqualityConstraint();

  /**
   * Destructor.
   */
  ~G06EqualityConstraint() override;

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
 * G06 constrained test problem.
 *
 * * Number of parameters: 2
 * * Number of inequality constraints: 2
 * * Number of equality constraints: 0
 * * Domain: \f$\bar{\vec{x}} \in [13, 100] \times [0, 100]\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (2819/200, 0.8429608)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -6961.814\f$
 */
class G06 : public ConstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  G06();

  /**
   * Destructor.
   */
  ~G06() override;

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
  G06Objective f;
  /// inequality constraint function
  G06InequalityConstraint g;
  /// equality constraint function
  G06EqualityConstraint h;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G06_HPP */

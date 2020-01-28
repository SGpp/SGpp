// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G05_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G05_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * G05 objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * 3 \bar{x}_1 + \left(\frac{\bar{x}_1}{100}\right)^3 +
 * \frac{2}{3} \left(\frac{\bar{x}_2}{100}\right)^3\f]
 */
class G05Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  G05Objective();

  /**
   * Destructor.
   */
  ~G05Objective() override;

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
 * G05 inequality constraint function.
 *
 * Definition:
 * \f[\bar{\vec{g}}(\bar{\vec{x}}) :=
 * \begin{pmatrix}
 *     -\bar{x}_4 + \bar{x}_3 - 11/20\\
 *     -\bar{x}_3 + \bar{x}_4 - 11/20
 * \end{pmatrix}\f]
 */
class G05InequalityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G05InequalityConstraint();

  /**
   * Destructor.
   */
  ~G05InequalityConstraint() override;

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
 * G05 equality constraint function.
 *
 * Definition:
 * \f[\bar{\vec{h}}(\bar{\vec{x}}) :=
 * \begin{pmatrix}
 *     1000 \sin(-\bar{x}_3 - 1/4) +
 *     1000 \sin(-\bar{x}_4 - 1/4) +
 *     894.8 - \bar{x}_1\\
 *     1000 \sin(\bar{x}_3 - 1/4) +
 *     1000 \sin(\bar{x}_3 - \bar{x}_4 - 1/4) +
 *     894.8 - \bar{x}_2\\
 *     1000 \sin(\bar{x}_4 - 1/4) +
 *     1000 \sin(\bar{x}_4 - \bar{x}_3 - 1/4) +
 *     1294.8
 * \end{pmatrix}\f]
 */
class G05EqualityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G05EqualityConstraint();

  /**
   * Destructor.
   */
  ~G05EqualityConstraint() override;

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
 * G05 constrained test problem.
 *
 * * Number of parameters: 4
 * * Number of inequality constraints: 2
 * * Number of equality constraints: 3
 * * Domain: \f$\bar{\vec{x}} \in [0, 1200]^2 \times [-0.55, 0.55]^2\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (679.9453, 1026.067, 0.1188764, -0.3962336)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   5126.498\f$
 */
class G05 : public ConstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  G05();

  /**
   * Destructor.
   */
  ~G05() override;

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
  G05Objective f;
  /// inequality constraint function
  G05InequalityConstraint g;
  /// equality constraint function
  G05EqualityConstraint h;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G05_HPP */

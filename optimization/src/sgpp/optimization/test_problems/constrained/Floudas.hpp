// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_FLOUDAS_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_FLOUDAS_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Floudas objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * -\bar{x}_1 - \bar{x}_2\f]
 */
class FloudasObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  FloudasObjective();

  /**
   * Destructor.
   */
  ~FloudasObjective() override;

  /**
   * @param x     point \f$\vec{x} \in [0, 1]^d\f$
   * @return      \f$f(\vec{x})\f$
   */
  double evalUndisplaced(const base::DataVector& x) override;

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<ScalarFunction>& clone) const override;
};

/**
 * Floudas inequality constraint function.
 *
 * Definition:
 * \f[\vec{\bar{g}}(\bar{\vec{x}}) :=
 * \begin{pmatrix}
 *     \bar{x}_2 - 2 \bar{x}_1^4 + 8 \bar{x}_1^3 - 8 \bar{x}_1^2 - 2\\
 *     \bar{x}_2 - 4 \bar{x}_1^4 + 32 \bar{x}_1^3 - 88 \bar{x}_1^2 +
 *     96 \bar{x}_1 - 36
 * \end{pmatrix}\f]
 */
class FloudasInequalityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  FloudasInequalityConstraint();

  /**
   * Destructor.
   */
  ~FloudasInequalityConstraint() override;

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
 * Floudas equality constraint function.
 *
 * Definition: empty, i.e., no constraint
 */
class FloudasEqualityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  FloudasEqualityConstraint();

  /**
   * Destructor.
   */
  ~FloudasEqualityConstraint() override;

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
 * Floudas constrained test problem.
 *
 * * Number of parameters: 2
 * * Number of inequality constraints: 2
 * * Number of equality constraints: 0
 * * Domain: \f$\bar{\vec{x}} \in [0, 3] \times [0, 4]\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (2.329520, 3.178493)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -5.508013\f$
 */
class Floudas : public ConstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  Floudas();

  /**
   * Destructor.
   */
  ~Floudas() override;

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
  FloudasObjective f;
  /// inequality constraint function
  FloudasInequalityConstraint g;
  /// equality constraint function
  FloudasEqualityConstraint h;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_FLOUDAS_HPP */

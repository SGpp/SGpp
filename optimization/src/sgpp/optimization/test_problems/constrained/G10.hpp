// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G10_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G10_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace SGPP {
namespace optimization {
namespace test_problems {

/**
 * G10 objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * \bar{x}_1 + \bar{x}_2 + \bar{x}_3\f]
 */
class G10Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  G10Objective();

  /**
   * Destructor.
   */
  virtual ~G10Objective() override;

  /**
   * @param x     point \f$\vec{x} \in [0, 1]^d\f$
   * @return      \f$f(\vec{x})\f$
   */
  virtual float_t evalUndisplaced(const base::DataVector& x)
  override;

  /**
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(std::unique_ptr<ScalarFunction>& clone)
  const override;
};

/**
 * G10 inequality constraint function.
 *
 * Definition:
 * \f[\bar{\vec{g}}(\bar{\vec{x}}) :=
 * \begin{pmatrix}
 *     -1 + (\bar{x}_4 + \bar{x}_6)/400\\
 *     -1 + (\bar{x}_5 + \bar{x}_7 - \bar{x}_4)/400\\
 *     -1 + (\bar{x}_8 - \bar{x}_5)/100\\
 *     -\bar{x}_1 \bar{x}_6 + 833.33252 \bar{x}_4 + 100 \bar{x}_1 -
 *     83333.333\\
 *     -\bar{x}_2 \bar{x}_7 + 1250 \bar{x}_5 + \bar{x}_2 \bar{x}_4 -
 *     1250 \bar{x}_4\\
 *     -\bar{x}_3 \bar{x}_8 + 1250000 + \bar{x}_3 \bar{x}_5 -
 *     2500 \bar{x}_5
 * \end{pmatrix}\f]
 */
class G10InequalityConstraint :
  public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G10InequalityConstraint();

  /**
   * Destructor.
   */
  virtual ~G10InequalityConstraint() override;

  /**
   * @param       x       point \f$\vec{x} \in \mathbb{R}^d\f$
   * @param[out]  value   \f$\vec{f}(\vec{x})\f$
   */
  virtual void evalUndisplaced(const base::DataVector& x,
                               base::DataVector& value) override;

  /**
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(std::unique_ptr<VectorFunction>& clone)
  const override;
};

/**
 * G10 equality constraint function.
 *
 * Definition: empty, i.e., no constraint
 */
class G10EqualityConstraint :
  public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G10EqualityConstraint();

  /**
   * Destructor.
   */
  virtual ~G10EqualityConstraint() override;

  /**
   * @param       x       point \f$\vec{x} \in \mathbb{R}^d\f$
   * @param[out]  value   \f$\vec{f}(\vec{x})\f$
   */
  virtual void evalUndisplaced(const base::DataVector& x,
                               base::DataVector& value) override;

  /**
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(std::unique_ptr<VectorFunction>& clone)
  const override;
};

/**
 * G10 constrained test problem.
 *
 * * Number of parameters: 8
 * * Number of inequality constraints: 6
 * * Number of equality constraints: 0
 * * Domain: \f$\bar{\vec{x}} \in
 *   [10^2, 10^4] \times [10^3, 10^4]^2 \times [10, 10^3]^5\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (579.3167, 1359.943, 5110.071, 182.0174,
 *   295.5985, 217.9799, 286.4162, 395.5979)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   7049.331\f$
 */
class G10 : public ConstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  G10();

  /**
   * Destructor.
   */
  virtual ~G10() override;

  /**
   * @return  objective function of the test problem
   */
  virtual TestScalarFunction& getObjectiveFunction() override;

  /**
   * @return  inequality function of the test problem
   */
  virtual TestVectorFunction& getInequalityConstraintFunction() override;

  /**
   * @return  equality constraint function of the test problem
   */
  virtual TestVectorFunction& getEqualityConstraintFunction() override;

  /**
   * @param[out] x minimal point
   *               \f$\vec{x}_\opt \in [0, 1]^d\f$
   * @return       minimal function value
   *               \f$f(\vec{x}_\opt)\f$
   */
  virtual float_t getOptimalPointUndisplaced(base::DataVector& x)
  override;

 protected:
  /// objective function
  G10Objective f;
  /// inequality constraint function
  G10InequalityConstraint g;
  /// equality constraint function
  G10EqualityConstraint h;
};

}
}
}

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G10_HPP */

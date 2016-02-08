// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_SIMIONESCU_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_SIMIONESCU_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace SGPP {
namespace optimization {
namespace test_problems {

/**
 * Simionescu objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * \frac{1}{10} \bar{x}_1 \bar{x}_2\f]
 */
class SimionescuObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  SimionescuObjective();

  /**
   * Destructor.
   */
  virtual ~SimionescuObjective() override;

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
 * Simionescu inequality constraint function.
 *
 * Definition:
 * \f[\bar{g}(\bar{\vec{x}}) :=
 * \bar{x}_1^2 + \bar{x}_2^2 - \left(1 + \frac{1}{5}
 * \cos\!\Big(8 \arctan(\bar{x}_1/\bar{x}_2)\Big)\right)\f]
 */
class SimionescuInequalityConstraint :
  public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  SimionescuInequalityConstraint();

  /**
   * Destructor.
   */
  virtual ~SimionescuInequalityConstraint() override;

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
 * Simionescu equality constraint function.
 *
 * Definition: empty, i.e., no constraint
 */
class SimionescuEqualityConstraint :
  public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  SimionescuEqualityConstraint();

  /**
   * Destructor.
   */
  virtual ~SimionescuEqualityConstraint() override;

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
 * Simionescu constrained test problem.
 *
 * * Number of parameters: 2
 * * Number of inequality constraints: 1
 * * Number of equality constraints: 0
 * * Domain: \f$\bar{\vec{x}} \in [-5/4, 5/4]^2\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   3\sqrt{2}/5 \cdot (\pm 1, \mp 1)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -9/125\f$
 */
class Simionescu : public ConstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  Simionescu();

  /**
   * Destructor.
   */
  virtual ~Simionescu() override;

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
  SimionescuObjective f;
  /// inequality constraint function
  SimionescuInequalityConstraint g;
  /// equality constraint function
  SimionescuEqualityConstraint h;
};

}
}
}

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_SIMIONESCU_HPP */

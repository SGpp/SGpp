// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G12_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G12_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace SGPP {
namespace optimization {
namespace test_problems {

/**
 * G12 objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * \norm{\bar{\vec{x}} - \vec{5}}_2^2/100 - 1\f]
 */
class G12Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  G12Objective();

  /**
   * Destructor.
   */
  virtual ~G12Objective() override;

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
 * G12 inequality constraint function.
 *
 * Definition:
 * \f[\bar{g}(\bar{\vec{x}}) :=
 * \min_{y_1,y_2,y_3 = 1,\dotsc,9}
 * \norm{\bar{\vec{x}} - \vec{y}}_2^2 - 1/16\f]
 */
class G12InequalityConstraint :
  public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G12InequalityConstraint();

  /**
   * Destructor.
   */
  virtual ~G12InequalityConstraint() override;

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
 * G12 equality constraint function.
 *
 * Definition: empty, i.e., no constraint
 */
class G12EqualityConstraint :
  public TestVectorFunction {
 public:
  /**
   * Constructor.
   */
  G12EqualityConstraint();

  /**
   * Destructor.
   */
  virtual ~G12EqualityConstraint() override;

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
 * G12 constrained test problem.
 *
 * * Number of parameters: 3
 * * Number of inequality constraints: 1
 * * Number of equality constraints: 0
 * * Domain: \f$\bar{\vec{x}} \in [0, 10]^3\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   (5, 5, 5)\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -1\f$
 */
class G12 : public ConstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  G12();

  /**
   * Destructor.
   */
  virtual ~G12() override;

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
  G12Objective f;
  /// inequality constraint function
  G12InequalityConstraint g;
  /// equality constraint function
  G12EqualityConstraint h;
};

}
}
}

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G12_HPP */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G03_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G03_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * G03 objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) :=
 * -\prod_{t=1}^d \bar{x}_t\f]
 */
class G03Objective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   *
   * @param d   number of parameters
   */
  explicit G03Objective(size_t d);

  /**
   * Destructor.
   */
  ~G03Objective() override;

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
 * G03 inequality constraint function.
 *
 * Definition: empty, i.e., no constraint
 */
class G03InequalityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   *
   * @param d   number of parameters
   */
  explicit G03InequalityConstraint(size_t d);

  /**
   * Destructor.
   */
  ~G03InequalityConstraint() override;

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
 * G03 equality constraint function.
 *
 * Definition:
 * \f[\bar{h}(\bar{\vec{x}}) :=
 * \norm{\bar{\vec{x}}}_2^2 - 1\f]
 */
class G03EqualityConstraint : public TestVectorFunction {
 public:
  /**
   * Constructor.
   *
   * @param d   number of parameters
   */
  explicit G03EqualityConstraint(size_t d);

  /**
   * Destructor.
   */
  ~G03EqualityConstraint() override;

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
 * G03 constrained test problem.
 *
 * * Number of parameters: \f$d\f$
 * * Number of inequality constraints: 0
 * * Number of equality constraints: 1
 * * Domain: \f$\bar{\vec{x}} \in [0, 1]^d\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   d^{-1/2} \cdot \vec{1}\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -d^{-d/2}\f$
 */
class G03 : public ConstrainedTestProblem {
 public:
  /**
   * Constructor.
   *
   * @param d   number of parameters
   */
  explicit G03(size_t d);

  /**
   * Destructor.
   */
  ~G03() override;

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
  G03Objective f;
  /// inequality constraint function
  G03InequalityConstraint g;
  /// equality constraint function
  G03EqualityConstraint h;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_CONSTRAINED_G03_HPP */

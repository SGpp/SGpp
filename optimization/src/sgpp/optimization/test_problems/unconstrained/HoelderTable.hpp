// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_HOELDERTABLE_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_HOELDERTABLE_HPP

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

namespace SGPP {
namespace optimization {
namespace test_problems {

/**
 * Hoelder table objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) := -\left|\sin \bar{x}_1 \cos \bar{x}_2
 * \exp\!\left(\left|1 -
 * \frac{\norm{\bar{\vec{x}}}_2}{\pi}\right|\right)\right|\f],
 */
class HoelderTableObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   */
  HoelderTableObjective();

  /**
   * Destructor.
   */
  virtual ~HoelderTableObjective() override;

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
 * Hoelder table unconstrained test problem.
 *
 * * Number of parameters: 2
 * * Domain: \f$\bar{\vec{x}} \in [-10, 10]^2\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} \in
 *   \{(8.055023, \pm 9.664590), (-8.055023, \pm 9.664590)\}\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   -19.20850\f$
 *
 * The displacement is restricted because the minimal points lie near
 * the corners of \f$[0, 1]^2\f$.
 */
class HoelderTable : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   */
  HoelderTable();

  /**
   * Destructor.
   */
  virtual ~HoelderTable() override;

  /**
   * @return  objective function of the test problem
   */
  virtual TestScalarFunction& getObjectiveFunction() override;

  /**
   * @param[out] x minimal point
   *               \f$\vec{x}_\opt \in [0, 1]^d\f$
   * @return       minimal function value
   *               \f$f(\vec{x}_\opt)\f$
   */
  virtual float_t getOptimalPointUndisplaced(base::DataVector& x)
  override;

 protected:
  /**
   * Checks if \f$\vec{d} \in [-0.005, 0.005] \times [-0.01, 0.01]\f$.
   *
   * @return whether the current displacement is feasible
   */
  virtual bool isDisplacementFeasible() override;

  /// objective function
  HoelderTableObjective f;
};

}
}
}

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_HOELDERTABLE_HPP */

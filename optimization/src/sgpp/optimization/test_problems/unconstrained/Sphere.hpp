// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_SPHERE_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_SPHERE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>

namespace SGPP {
namespace optimization {
namespace test_problems {

/**
 * Sphere objective function.
 *
 * Definition:
 * \f[\bar{f}(\bar{\vec{x}}) := \norm{\bar{\vec{x}}}_2^2\f]
 */
class SphereObjective : public TestScalarFunction {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  SphereObjective(size_t d);

  /**
   * Destructor.
   */
  virtual ~SphereObjective() override;

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
 * Sphere unconstrained test problem.
 *
 * * Number of parameters: \f$d\f$
 * * Domain: \f$\bar{\vec{x}} \in [-1, 9]^d\f$
 * * Optimal point: \f$\bar{\vec{x}}_{\text{opt}} =
 *   \vec{0}\f$
 * * Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =
 *   0\f$
 */
class Sphere : public UnconstrainedTestProblem {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  Sphere(size_t d);

  /**
   * Destructor.
   */
  virtual ~Sphere() override;

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
  /// objective function
  SphereObjective f;
};

}
}
}

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_SPHERE_HPP */

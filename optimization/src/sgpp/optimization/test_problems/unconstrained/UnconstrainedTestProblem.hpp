// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_UNCONSTRAINEDTESTPROBLEM_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_UNCONSTRAINEDTESTPROBLEM_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/test_problems/TestScalarFunction.hpp>

#include <cstddef>

namespace SGPP {
namespace optimization {
namespace test_problems {

/**
 * Base class for analytical, unconstrained test problems.
 * This class essentially manages an objective function,
 * generates a Gaussian displacement vector, and
 * contains the location of the optimal point.
 */
class UnconstrainedTestProblem {
 public:
  /// default standard deviation for the displacement vector
  static constexpr float_t DEFAULT_STANDARD_DEVIATION = 0.01;

  /**
   * Constructor.
   * The displacement is set to all zeros, so to displace the function
   * call generateDisplacement() afterwards.
   *
   * @param d     dimension of the domain
   */
  explicit UnconstrainedTestProblem(size_t d);

  /**
   * Destructor.
   */
  virtual ~UnconstrainedTestProblem();

  /**
   * @return  objective function of the test problem
   */
  virtual TestScalarFunction& getObjectiveFunction() = 0;

  /**
   * Returns the minimal point of the displaced function.
   *
   * @param[out] x reverse displaced minimal point
   *               \f$\vec{x}_\opt - \vec{d}\f$
   * @return       minimal function value
   *               \f$f(\vec{x}_\opt)\f$
   */
  float_t getOptimalPoint(base::DataVector& x);

  /**
   * Pure virtual method returning the minimal point
   * (or one of the minimal points, if there are multiple of them)
   * of the test function.
   *
   * @param[out] x    minimal point
   *                  \f$\vec{x}_\opt\f$
   * @return          minimal function value
   *                  \f$f(\vec{x}_\opt)\f$
   */
  virtual float_t getOptimalPointUndisplaced(base::DataVector& x) = 0;

  /**
   * Generate normally distributed pseudorandom displacement with
   * default standard deviation.
   * This method also sets the new displacement in the
   * objective function.
   */
  void generateDisplacement();

  /**
   * Generate normally distributed pseudorandom displacement.
   * This method also sets the new displacement in the
   * objective function.
   *
   * @param stdDev standard deviation of the displacement coordinates
   */
  void generateDisplacement(float_t stdDev);

  /**
   * @return                currently used displacement
   */
  const base::DataVector& getDisplacement() const;

  /**
   * Sets the displacement vector.
   * This method also sets the new displacement in the
   * objective function.
   *
   * @param displacement    currently used displacement
   */
  void setDisplacement(const base::DataVector& displacement);

 protected:
  /// number of parameters
  size_t d;
  /// displacement vector
  base::DataVector displacement;

  /**
   * Checks if the current displacement is ok for the specific
   * objective function.
   * An infeasible displacement would be one for which one of the
   * following happens:
   * 1. The optimal point is displaced out of \f$[0, 1]^d\f$.
   * 2. Singularities, where the function is not well-defined,
   *    are displaced into \f$[0, 1]^d\f$.
   *
   * This function is called by generateDisplacement
   * to check the feasibility of the new displacement.
   * Only 1. will be checked by default.
   * For additional constraints on the displacement, override this
   * function.
   *
   * @return whether the current displacement is feasible
   */
  virtual bool isDisplacementFeasible();
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace SGPP

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_UNCONSTRAINED_UNCONSTRAINEDTESTPROBLEM_HPP */

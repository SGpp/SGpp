// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_TESTSCALARFUNCTION_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_TESTSCALARFUNCTION_HPP

#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Base class for analytical objective function examples
 * ("test functions").
 * This class differs from ScalarFunction in the way that
 * evaluation point \f$\vec{x}\f$ are pseudorandomly displaced
 * when calling the eval() function.
 * Evaluting the undisplaced function must be implemented in
 * evalUndisplaced.
 * The displaced function is
 * \f$\vec{x} \mapsto f(\vec{x} + \vec{d})\f$
 * for a vector \f$\vec{d}\f$ ("displacement").
 *
 * Taking the average of results of multiple runs with different
 * displacements makes results more robust and significant.
 */
class TestScalarFunction : public base::ScalarFunction {
 public:
  /**
   * Constructor.
   * The displacement is set to all zeros, so to displace the function
   * call generateDisplacement() afterwards.
   *
   * @param d     dimension of the domain
   */
  explicit TestScalarFunction(size_t d);

  /**
   * Destructor.
   */
  ~TestScalarFunction() override;

  /**
   * Evaluate displaced function.
   *
   * @param x point \f$\vec{x} \in \mathbb{R}^d\f$
   * @return  \f$f(\vec{x} + \vec{d})\f$
   *          with displacement \f$\vec{d}\f$
   */
  double eval(const base::DataVector& x) override;

  /**
   * Pure virtual method for evaluating the undisplaced function.
   *
   * @param x     point \f$\vec{x} \in \mathbb{R}^d\f$
   * @return      \f$f(\vec{x})\f$
   */
  virtual double evalUndisplaced(const base::DataVector& x) = 0;

  /**
   * @return                currently used displacement
   */
  const base::DataVector& getDisplacement() const;

  /**
   * @param displacement    currently used displacement
   */
  void setDisplacement(const base::DataVector& displacement);

 protected:
  /// displacement vector
  base::DataVector displacement;
  /// temporary vector for displacing
  base::DataVector xTmp;
};
}  // namespace test_problems
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_TESTSCALARFUNCTION_HPP */

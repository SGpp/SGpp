// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TEST_PROBLEMS_TESTVECTORFUNCTION_HPP
#define SGPP_OPTIMIZATION_TEST_PROBLEMS_TESTVECTORFUNCTION_HPP

#include <sgpp/base/function/vector/VectorFunction.hpp>

#include <sgpp/globaldef.hpp>
#include <cstddef>

namespace sgpp {
namespace optimization {
namespace test_problems {

/**
 * Base class for analytical constraint function examples
 * ("test functions").
 * This class differs from VectorFunction in the way that
 * evaluation point \f$\vec{x}\f$ are pseudorandomly displaced
 * when calling the eval() function.
 * Evaluting the undisplaced function must be implemented in
 * evalUndisplaced.
 * The displaced function is
 * \f$\vec{x} \mapsto \vec{f}(\vec{x} + \vec{d})\f$
 * for a vector \f$\vec{d}\f$ ("displacement").
 */
class TestVectorFunction : public base::VectorFunction {
 public:
  /**
   * Constructor.
   * The displacement is set to all zeros, so to displace the function
   * call generateDisplacement() afterwards.
   *
   * @param d     dimension of the domain
   * @param m     number of components
   */
  TestVectorFunction(size_t d, size_t m);

  /**
   * Destructor.
   */
  ~TestVectorFunction() override;

  /**
   * Evaluate displaced function.
   *
   * @param       x       point \f$\vec{x} \in \mathbb{R}^d\f$
   * @param[out]  value   \f$\vec{f}(\vec{x} + \vec{d})\f$
   *                      with displacement \f$\vec{d}\f$
   */
  void eval(const base::DataVector& x, base::DataVector& value) override;

  /**
   * Pure virtual method for evaluating the undisplaced function.
   *
   * @param       x       point \f$\vec{x} \in \mathbb{R}^d\f$
   * @param[out]  value   \f$\vec{f}(\vec{x})\f$
   */
  virtual void evalUndisplaced(const base::DataVector& x, base::DataVector& value) = 0;

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

#endif /* SGPP_OPTIMIZATION_TEST_PROBLEMS_TESTVECTORFUNCTION_HPP */

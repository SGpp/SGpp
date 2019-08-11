// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <memory>

namespace sgpp {
namespace base {

/**
 * Abstract base class for scalar-valued functions
 * \f$f\colon [0, 1]^d \to \mathbb{R}\f$
 * together with their gradients
 * \f$\nabla f\colon [0, 1]^d \to \mathbb{R}^d\f$
 * (e.g., gradients of objective functions in optimization).
 */
class ScalarFunctionGradient {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  explicit ScalarFunctionGradient(size_t d) : d(d) {}

  /**
   * Destructor.
   */
  virtual ~ScalarFunctionGradient() {}

  /**
   * Pure virtual method for calculating
   * \f$f(\vec{x})\f$ together with \f$\nabla f(\vec{x})\f$.
   *
   * @param      x        evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @param[out] gradient gradient
   *                      \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
   * @return              \f$f(\vec{x})\f$
   */
  virtual double eval(const DataVector& x, DataVector& gradient) = 0;

  /**
   * Convenience method for calculating \f$f(\vec{x})\f$ together with
   * \f$\nabla f(\vec{x})\f$ for multiple \f$\vec{x}\f$.
   *
   * @param      x        matrix \f$\vec{x} \in [0, 1]^{N \times d}\f$
   *                      of evaluation points (row-wise)
   * @param[out] value    vector of size \f$N\f$, where the \f$k\f$-th
   *                      entry is \f$f(\vec{x}_k)\f$
   *                      (where \f$\vec{x}_k\f$ is the \f$k\f$-th row
   *                      of \f$x\f$)
   * @param[out] gradient matrix of size \f$N \times d\f$
   *                      where the \f$k\f$-th row is
   *                      \f$\nabla f(\vec{x}_k)\f$
   */
  virtual void eval(const DataMatrix& x, DataVector& value,
                    DataMatrix& gradient) {
    const size_t N = x.getNrows();
    DataVector xk(d);
    DataVector yk(d);
    value.resize(N);
    gradient.resize(N, d);

    for (size_t k = 0; k < N; k++) {
      x.getRow(k, xk);
      value[k] = eval(xk, yk);
      gradient.setRow(k, yk);
    }
  }


  /**
   * @return dimension \f$d\f$ of the domain
   */
  size_t getNumberOfParameters() const { return d; }

  /**
   * Pure virtual method for cloning the gradient.
   * It should generate a pointer to the cloned object and
   * it's used for parallel computations
   * (the eval() method might not be thread-safe).
   *
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const = 0;

 protected:
  /// dimension of the domain
  size_t d;
};
}  // namespace base
}  // namespace sgpp

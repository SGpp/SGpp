// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <memory>

namespace sgpp {
namespace base {

/**
 * Abstract base class for scalar-valued functions
 * \f$f\colon [0, 1]^d \to \mathbb{R}\f$
 * (e.g., objective functions in optimization).
 */
class ScalarFunction {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  explicit ScalarFunction(size_t d) : d(d) {}

  /**
   * Destructor.
   */
  virtual ~ScalarFunction() {}

  /**
   * Pure virtual method for calculating \f$f(\vec{x})\f$.
   *
   * @param x     evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @return      \f$f(\vec{x})\f$
   */
  virtual double eval(const DataVector& x) = 0;

  /**
   * Convenience method for calculating \f$f(\vec{x})\f$ for multiple \f$\vec{x}\f$.
   *
   * @param      x      matrix \f$\vec{x} \in [0, 1]^{N \times d}\f$
   *                    of evaluation points (row-wise)
   * @param[out] value  \f$(f(\vec{x}_k))_k\f$
   *                    where \f$\vec{x}_k\f$ is the \f$k\f$-th row of \f$x\f$
   */
  virtual void eval(const DataMatrix& x, DataVector& value) {
    const size_t N = x.getNrows();
    DataVector xk(d);
    value.resize(N);

    for (size_t k = 0; k < N; k++) {
      x.getRow(k, xk);
      value[k] = eval(xk);
    }
  }

  /**
   * @return dimension \f$d\f$ of the domain
   */
  size_t getNumberOfParameters() const { return d; }

  /**
   * Pure virtual method for cloning the function.
   * It should generate a pointer to the cloned object and
   * it's used for parallel computations
   * (the eval() method might not be thread-safe).
   *
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(std::unique_ptr<ScalarFunction>& clone) const = 0;

 protected:
  /// dimension of the domain
  size_t d;
};
}  // namespace base
}  // namespace sgpp

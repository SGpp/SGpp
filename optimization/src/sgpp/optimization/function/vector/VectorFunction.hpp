// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_VECTOR_VECTORFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_VECTOR_VECTORFUNCTION_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <cstddef>
#include <memory>

namespace sgpp {
namespace optimization {

/**
 * Abstract base class for vector-valued functions
 * \f$g\colon [0, 1]^d \to \mathbb{R}^m\f$
 * (e.g., equality/inequality constraints
 * \f$g(\vec{x}) \le \vec{0}\f$ or \f$g(\vec{x}) = \vec{0}\f$
 * in optimization).
 */
class VectorFunction {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   * @param m     number of components
   */
  VectorFunction(size_t d, size_t m) : d(d), m(m) {}

  /**
   * Destructor.
   */
  virtual ~VectorFunction() {}

  /**
   * Pure virtual method for calculating \f$g(\vec{x})\f$.
   *
   * @param[in]  x      evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @param[out] value  \f$g(\vec{x})\f$
   */
  virtual void eval(const base::DataVector& x, base::DataVector& value) = 0;

  /**
   * Convenience method for calculating \f$g(\vec{x})\f$ for multiple \f$\vec{x}\f$.
   *
   * @param      x      matrix \f$\vec{x} \in [0, 1]^{N \times d}\f$
   *                    of evaluation points (row-wise)
   * @param[out] value  matrix of size \f$N \times m\f$
   *                    where the \f$k\f$-th row is \f$g(\vec{x}_k)\f$
   *                    (where \f$\vec{x}_k\f$ is the \f$k\f$-th row of \f$x\f$)
   */
  virtual void eval(const base::DataMatrix& x, base::DataMatrix& value) {
    const size_t N = x.getNrows();
    base::DataVector xk(d);
    base::DataVector yk(m);
    value.resize(N, m);

    for (size_t k = 0; k < N; k++) {
      x.getRow(k, xk);
      eval(xk, yk);
      value.setRow(k, yk);
    }
  }

  /**
   * @return dimension \f$d\f$ of the domain
   */
  size_t getNumberOfParameters() const { return d; }

  /**
   * @return number \f$m\f$ of components
   */
  size_t getNumberOfComponents() const { return m; }

  /**
   * Pure virtual method for cloning the function.
   * It should generate a pointer to the cloned object and
   * it's used for parallel computations
   * (the eval() method might not be thread-safe).
   *
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(std::unique_ptr<VectorFunction>& clone) const = 0;

 protected:
  /// dimension of the domain
  size_t d;
  /// number of components
  size_t m;
};
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUNCTION_VECTOR_VECTORFUNCTION_HPP */

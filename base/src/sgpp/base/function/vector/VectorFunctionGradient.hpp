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
 * Abstract base class for vector-valued functions
 * \f$g\colon [0, 1]^d \to \mathbb{R}^m\f$
 * together with their Jacobians
 * \f$\nabla g\colon [0, 1]^d \to \mathbb{R}^{m \times d}\f$,
 * i.e. \f$(\nabla g)_{i,t} = \frac{\partial}{\partial x_t} g_i\f$
 * with \f$g = (g_i)_{i=1}^m\f$
 * (e.g., equality/inequality constraints
 * \f$g(\vec{x}) \le \vec{0}\f$ or \f$g(\vec{x}) = \vec{0}\f$
 * in optimization).
 */
class VectorFunctionGradient {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   * @param m     number of components
   */
  VectorFunctionGradient(size_t d, size_t m) : d(d), m(m) {}

  /**
   * Destructor.
   */
  virtual ~VectorFunctionGradient() {}

  /**
   * Pure virtual method for calculating \f$g(\vec{x})\f$
   * together with \f$\nabla g(\vec{x})\f$.
   *
   * @param[in]  x        evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @param[out] value    \f$g(\vec{x})\f$
   * @param[out] gradient gradient \f$\nabla g(\vec{x}) \in
   *                      \mathbb{R}^{m \times d}\f$
   */
  virtual void eval(const DataVector& x, DataVector& value,
                    DataMatrix& gradient) = 0;

  /**
   * @return dimension \f$d\f$ of the domain
   */
  size_t getNumberOfParameters() const { return d; }

  /**
   * @return number \f$m\f$ of components
   */
  size_t getNumberOfComponents() const { return m; }

  /**
   * Pure virtual method for cloning the gradient.
   * It should generate a pointer to the cloned object and
   * it's used for parallel computations
   * (the eval() method might not be thread-safe).
   *
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(std::unique_ptr<VectorFunctionGradient>& clone) const = 0;

 protected:
  /// dimension of the domain
  size_t d;
  /// number of components
  size_t m;
};
}  // namespace base
}  // namespace sgpp

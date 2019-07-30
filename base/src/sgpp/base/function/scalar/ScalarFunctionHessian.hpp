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
#include <vector>

namespace sgpp {
namespace base {

/**
 * Abstract base class for scalar-valued functions
 * \f$f\colon [0, 1]^d \to \mathbb{R}\f$
 * together with their gradients
 * \f$\nabla f\colon [0, 1]^d \to \mathbb{R}^d\f$
 * and Hessians
 * \f$H_f\colon [0, 1]^d \to \mathbb{R}^{d \times d}\f$
 * (e.g., Hessians of objective functions in optimization).
 */
class ScalarFunctionHessian {
 public:
  /**
   * Constructor.
   *
   * @param d     dimension of the domain
   */
  explicit ScalarFunctionHessian(size_t d) : d(d) {}

  /**
   * Destructor.
   */
  virtual ~ScalarFunctionHessian() {}

  /**
   * Pure virtual method for calculating \f$f(\vec{x})\f$ together with
   * \f$\nabla f(\vec{x})\f$ and
   * \f$H_f(\vec{x}) \in \mathbb{R}^{d \times d}\f$.
   *
   * @param      x        evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @param[out] gradient gradient
   *                      \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
   * @param[out] hessian  Hessian matrix
   *                      \f$H_f(\vec{x}) \in \mathbb{R}^{d \times d}\f$
   * @return              \f$f(\vec{x})\f$
   */
  virtual double eval(const base::DataVector& x, base::DataVector& gradient,
                      base::DataMatrix& hessian) = 0;

  /**
   * @return dimension \f$d\f$ of the domain
   */
  size_t getNumberOfParameters() const { return d; }

  /**
   * Pure virtual method for cloning the Hessian.
   * It should generate a pointer to the cloned object and
   * it's used for parallel computations
   * (the eval() method might not be thread-safe).
   *
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const = 0;

 protected:
  /// dimension of the domain
  size_t d;
};
}  // namespace base
}  // namespace sgpp

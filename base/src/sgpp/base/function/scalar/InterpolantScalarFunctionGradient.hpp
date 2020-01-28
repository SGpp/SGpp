// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradient.hpp>

#include <limits>

namespace sgpp {
namespace base {

/**
 * Sparse grid interpolant gradient of a scalar-valued function.
 *
 * @see InterpolantScalarFunction
 */
class InterpolantScalarFunctionGradient : public ScalarFunctionGradient {
 public:
  /**
   * Constructor.
   * Do not destruct the grid before the
   * InterpolantScalarFunctionGradient object!
   *
   * @param grid  sparse grid
   * @param alpha coefficient vector
   */
  InterpolantScalarFunctionGradient(Grid& grid, const DataVector& alpha)
      : ScalarFunctionGradient(grid.getDimension()),
        grid(grid),
        opEvalGradient(op_factory::createOperationEvalGradientNaive(grid)),
        alpha(alpha) {}

  /**
   * Destructor.
   */
  ~InterpolantScalarFunctionGradient() override {}

  /**
   * Evaluation of the function and its gradient.
   *
   * @param      x        evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @param[out] gradient gradient
   *                      \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
   * @return              \f$f(\vec{x})\f$
   */
  inline double eval(const DataVector& x, DataVector& gradient) override {
    for (size_t t = 0; t < d; t++) {
      if ((x[t] < 0.0) || (x[t] > 1.0)) {
        return std::numeric_limits<double>::infinity();
      }
    }

    return opEvalGradient->evalGradient(alpha, x, gradient);
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const override {
    clone =
        std::unique_ptr<ScalarFunctionGradient>(new InterpolantScalarFunctionGradient(grid, alpha));
  }

  /**
   * @return coefficient vector
   */
  const DataVector& getAlpha() const { return alpha; }

  /**
   * @param alpha coefficient vector
   */
  void setAlpha(const DataVector& alpha) { this->alpha = alpha; }

 protected:
  /// sparse grid
  Grid& grid;
  /// pointer to evaluation operation
  std::unique_ptr<OperationEvalGradient> opEvalGradient;
  /// coefficient vector
  DataVector alpha;
};
}  // namespace base
}  // namespace sgpp

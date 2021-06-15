// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>

#include <cstring>
#include <memory>

namespace sgpp {
namespace optimization {

/**
 * Sparse grid interpolant of a scalar-valued function.
 * Special case for active subspace response surfaces which is defined on mathbb{R}^d and not only
 * [0,1]^d
 *
 * More generally, the function can be any linear combination
 * $f colon mathbb{R}^d to mathbb{R}$,
 * $f(vec{x}) = sum_{k=1}^N alpha_k varphi_k(vec{x})$
 * of the basis functions
 * $varphi_k = varphi_{vec{ell}_k,vec{i}_k}$
 * of a sparse grid with grid points
 * $vec{x}_k = vec{x}_{vec{ell}_k,vec{i}_k}$.
 * But most often, the function (e.g., its coefficients) is constructed
 * as an interpolant at the grid points for some function values.
 */
class ASInterpolantScalarFunction : public base::ScalarFunction {
 public:
  /**
   * Constructor.
   * Do not destruct the grid before the
   * ASInterpolantScalarFunction object!
   *
   * @param grid  sparse grid
   * @param alpha coefficient vector
   */
  ASInterpolantScalarFunction(base::Grid& grid, const base::DataVector& alpha)
      : base::ScalarFunction(grid.getDimension()),
        grid(grid),
        opEval(op_factory::createOperationEvalNaive(grid)),
        alpha(alpha) {}

  /**
   * Destructor.
   */
  ~ASInterpolantScalarFunction() override {}

  /**
   * Evaluation of the function.
   *
   * @param x     evaluation point $vec{x} in mathbb{R}^d$
   * @return      $f(vec{x})$
   */
  inline double eval(const base::DataVector& x) override { return opEval->eval(alpha, x); }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<base::ScalarFunction>& clone) const override {
    clone = std::unique_ptr<base::ScalarFunction>(new ASInterpolantScalarFunction(grid, alpha));
  }

  /**
   * @return coefficient vector
   */
  const base::DataVector& getAlpha() const { return alpha; }

  /**
   * @param alpha coefficient vector
   */
  void setAlpha(const base::DataVector& alpha) { this->alpha = alpha; }

  /**
   * @return number of grid points
   */
  size_t getSize() { return this->grid.getSize(); }

 protected:
  /// sparse grid
  base::Grid& grid;
  /// pointer to evaluation operation
  std::unique_ptr<base::OperationEval> opEval;
  /// coefficient vector
  base::DataVector alpha;
};
}  // namespace optimization
}  // namespace sgpp

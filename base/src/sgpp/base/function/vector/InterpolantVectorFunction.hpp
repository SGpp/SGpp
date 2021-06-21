// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>

#include <cstddef>
#include <limits>

namespace sgpp {
namespace base {

/**
 * Sparse grid interpolant of a vector-valued function.
 *
 * More generally, the function can be any linear combination
 * \f$g\colon [0, 1]^d \to \mathbb{R}^m\f$,
 * \f$g_j(\vec{x}) = \sum_{k=1}^N \alpha_{k,j} \varphi_k(\vec{x})\f$
 * of the basis functions
 * \f$\varphi_k = \varphi_{\vec{\ell}_k,\vec{i}_k}\f$
 * of a sparse grid with grid points
 * \f$\vec{x}_k = \vec{x}_{\vec{\ell}_k,\vec{i}_k}\f$.
 * But most often, the function (e.g., its coefficients) is constructed
 * as an interpolant at the grid points for some function values.
 */
class InterpolantVectorFunction : public VectorFunction {
 public:
  /**
   * Constructor.
   * Do not destruct the grid before the
   * InterpolantVectorFunction object!
   *
   * @param grid  sparse grid
   * @param alpha coefficient matrix
   *              (j-th column contains hierarchical surplusses
   *              \f$\alpha_{\cdot,j}\f$ of \f$g_j\f$)
   */
  InterpolantVectorFunction(Grid& grid, const DataMatrix& alpha)
      : VectorFunction(grid.getDimension(), alpha.getNcols()),
        grid(grid),
        opEval(op_factory::createOperationEvalNaive(grid)),
        alpha(alpha) {}

  /**
   * Destructor.
   */
  ~InterpolantVectorFunction() override {}

  /**
   * Evaluation of the function.
   *
   * @param[in]  x      evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @param[out] value  \f$g(\vec{x})\f$
   */
  inline void eval(const DataVector& x, DataVector& value) override {
    for (size_t t = 0; t < d; t++) {
      if ((x[t] < 0.0) || (x[t] > 1.0)) {
        value.setAll(std::numeric_limits<double>::infinity());
        return;
      }
    }
    opEval->eval(alpha, x, value);
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<VectorFunction>& clone) const override {
    clone = std::unique_ptr<VectorFunction>(new InterpolantVectorFunction(grid, alpha));
  }

  /**
   * @return coefficient matrix
   */
  const DataMatrix& getAlpha() const { return alpha; }

  /**
   * @param alpha coefficient matrix
   */
  void setAlpha(const DataMatrix& alpha) { this->alpha = alpha; }

 protected:
  /// sparse grid
  Grid& grid;
  /// pointer to evaluation operation
  std::unique_ptr<OperationEval> opEval;
  /// coefficient matrix
  DataMatrix alpha;
};
}  // namespace base
}  // namespace sgpp

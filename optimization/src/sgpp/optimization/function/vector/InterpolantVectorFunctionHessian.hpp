// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_VECTOR_INTERPOLANTVECTORFUNCTIONHESSIAN_HPP
#define SGPP_OPTIMIZATION_FUNCTION_VECTOR_INTERPOLANTVECTORFUNCTIONHESSIAN_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradient.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionHessian.hpp>

#include <cstddef>
#include <vector>

namespace sgpp {
namespace optimization {

/**
 * Sparse grid interpolant Hessian of a vector-valued function.
 *
 * @see InterpolantVectorFunction
 */
class InterpolantVectorFunctionHessian : public VectorFunctionHessian {
 public:
  /**
   * Constructor.
   * Do not destruct the grid before the
   * InterpolantVectorFunctionGradient object!
   *
   * @param grid  sparse grid
   * @param alpha coefficient matrix
   *              (j-th column contains hierarchical surplusses
   *              \f$\alpha_{\cdot,j}\f$ of \f$g_j\f$)
   */
  InterpolantVectorFunctionHessian(base::Grid& grid, const base::DataMatrix& alpha)
      : VectorFunctionHessian(grid.getDimension(), alpha.getNcols()),
        grid(grid),
        opEvalHessian(op_factory::createOperationEvalHessianNaive(grid)),
        alpha(alpha) {}

  /**
   * Destructor.
   */
  ~InterpolantVectorFunctionHessian() override {}

  /**
   * Evaluation of the function and its gradient.
   *
   * @param[in]  x        evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @param[out] value    \f$g(\vec{x})\f$
   * @param[out] gradient Jacobian \f$\nabla g(\vec{x}) \in
   *                      \mathbb{R}^{m \times d}\f$
   * @param[out] hessian  \f$m\f$-vector of Hessians
   *                      \f$\nabla^2 g_i(\vec{x}) \in
   *                      \mathbb{R}^{d \times d}\f$
   */
  inline void eval(const base::DataVector& x, base::DataVector& value,
                   base::DataMatrix& gradient,
                   std::vector<base::DataMatrix>& hessian) override {
    for (size_t t = 0; t < d; t++) {
      if ((x[t] < 0.0) || (x[t] > 1.0)) {
        value.setAll(INFINITY);
        return;
      }
    }

    opEvalHessian->evalHessian(alpha, x, value, gradient, hessian);
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<VectorFunctionHessian>& clone) const override {
    clone =
        std::unique_ptr<VectorFunctionHessian>(new InterpolantVectorFunctionHessian(grid, alpha));
  }

  /**
   * @return coefficient matrix
   */
  const base::DataMatrix& getAlpha() const { return alpha; }

  /**
   * @param alpha coefficient matrix
   */
  void setAlpha(const base::DataMatrix& alpha) { this->alpha = alpha; }

 protected:
  /// sparse grid
  base::Grid& grid;
  /// pointer to evaluation operation
  std::unique_ptr<base::OperationEvalHessian> opEvalHessian;
  /// coefficient matrix
  base::DataMatrix alpha;
};
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUNCTION_VECTOR_INTERPOLANTVECTORFUNCTIONHESSIAN_HPP */

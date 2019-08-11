// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/function/vector/VectorFunctionHessian.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>
namespace sgpp {
namespace base {

/**
 * Implementation of VectorFunctionHessian that
 * wraps a std::function object.
 */
class WrapperVectorFunctionHessian : public VectorFunctionHessian {
 public:
  typedef std::function<void(const DataVector&, DataVector&, DataMatrix&,
                             std::vector<DataMatrix>&)>
      FunctionHessianEvalType;

  /**
   * Constructor.
   *
   * @param d         dimension of the domain
   * @param m         number of components
   * @param fHessian  function gradient to be wrapped
   */
  WrapperVectorFunctionHessian(size_t d, size_t m, FunctionHessianEvalType fHessian)
      : VectorFunctionHessian(d, m), fHessian(fHessian) {}

  /**
   * Destructor.
   */
  ~WrapperVectorFunctionHessian() override {}

  /**
   * @param[in]  x        evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @param[out] value    \f$g(\vec{x})\f$
   * @param[out] gradient Jacobian \f$\nabla g(\vec{x}) \in
   *                      \mathbb{R}^{m \times d}\f$
   * @param[out] hessian  \f$m\f$-vector of Hessians
   *                      \f$\nabla^2 g_i(\vec{x}) \in
   *                      \mathbb{R}^{d \times d}\f$
   */
  inline void eval(const DataVector& x, DataVector& value, DataMatrix& gradient,
                   std::vector<DataMatrix>& hessian) override {
    fHessian(x, value, gradient, hessian);
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<VectorFunctionHessian>& clone) const override {
    clone =
        std::unique_ptr<VectorFunctionHessian>(new WrapperVectorFunctionHessian(d, m, fHessian));
  }

 protected:
  /// function Hessian to be wrapped
  FunctionHessianEvalType fHessian;
};
}  // namespace base
}  // namespace sgpp

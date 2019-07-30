// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/function/vector/VectorFunctionGradient.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <functional>
#include <memory>
namespace sgpp {
namespace base {

/**
 * Implementation of VectorFunctionGradient that
 * wraps a std::function object.
 */
class WrapperVectorFunctionGradient : public VectorFunctionGradient {
 public:
  typedef std::function<void(const base::DataVector&, base::DataVector&, base::DataMatrix&)>
      FunctionGradientEvalType;

  /**
   * Constructor.
   *
   * @param d         dimension of the domain
   * @param m         number of components
   * @param fGradient function gradient to be wrapped
   */
  WrapperVectorFunctionGradient(size_t d, size_t m, FunctionGradientEvalType fGradient)
      : VectorFunctionGradient(d, m), fGradient(fGradient) {}

  /**
   * Destructor.
   */
  ~WrapperVectorFunctionGradient() override {}

  /**
   * @param[in]  x        evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @param[out] value    \f$g(\vec{x})\f$
   * @param[out] gradient Jacobian \f$\nabla g(\vec{x}) \in
   *                      \mathbb{R}^{m \times d}\f$
   */
  inline void eval(const base::DataVector& x, base::DataVector& value,
                   base::DataMatrix& gradient) override {
    fGradient(x, value, gradient);
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<VectorFunctionGradient>& clone) const override {
    clone =
        std::unique_ptr<VectorFunctionGradient>(new WrapperVectorFunctionGradient(d, m, fGradient));
  }

 protected:
  /// function gradient to be wrapped
  FunctionGradientEvalType fGradient;
};
}  // namespace base
}  // namespace sgpp

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <functional>
#include <memory>
namespace sgpp {
namespace base {

/**
 * Implementation of ScalarFunctionHessian that
 * wraps a std::function object.
 */
class WrapperScalarFunctionHessian : public ScalarFunctionHessian {
 public:
  typedef std::function<double(const base::DataVector&, base::DataVector&, base::DataMatrix&)>
      FunctionHessianEvalType;

  /**
   * Constructor.
   *
   * @param d         dimension of the domain
   * @param fHessian  function gradient to be wrapped
   */
  WrapperScalarFunctionHessian(size_t d, FunctionHessianEvalType fHessian)
      : ScalarFunctionHessian(d), fHessian(fHessian) {}

  /**
   * Destructor.
   */
  ~WrapperScalarFunctionHessian() override {}

  /**
   * @param      x        evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @param[out] gradient gradient
   *                      \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
   * @param[out] hessian  Hessian matrix
   *                      \f$H_f(\vec{x}) \in \mathbb{R}^{d \times d}\f$
   * @return              \f$f(\vec{x})\f$
   */
  inline double eval(const base::DataVector& x, base::DataVector& gradient,
                     base::DataMatrix& hessian) override {
    return fHessian(x, gradient, hessian);
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const override {
    clone = std::unique_ptr<ScalarFunctionHessian>(new WrapperScalarFunctionHessian(d, fHessian));
  }

 protected:
  /// function Hessian to be wrapped
  FunctionHessianEvalType fHessian;
};
}  // namespace base
}  // namespace sgpp

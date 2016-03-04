// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_SCALAR_WRAPPERSCALARFUNCTIONGRADIENT_HPP
#define SGPP_OPTIMIZATION_FUNCTION_SCALAR_WRAPPERSCALARFUNCTIONGRADIENT_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>

#include <cstddef>
#include <memory>
#include <functional>

namespace sgpp {
namespace optimization {

/**
 * Implementation of ScalarFunctionGradient that
 * wraps a std::function object.
 */
class WrapperScalarFunctionGradient : public ScalarFunctionGradient {
 public:
  typedef std::function<double(const base::DataVector&, base::DataVector&)>
      FunctionGradientEvalType;

  /**
   * Constructor.
   *
   * @param d         dimension of the domain
   * @param fGradient function gradient to be wrapped
   */
  WrapperScalarFunctionGradient(size_t d, FunctionGradientEvalType fGradient)
      : ScalarFunctionGradient(d), fGradient(fGradient) {}

  /**
   * Destructor.
   */
  ~WrapperScalarFunctionGradient() override {}

  /**
   * @param      x        evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @param[out] gradient gradient
   *                      \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
   * @return              \f$f(\vec{x})\f$
   */
  inline double eval(const base::DataVector& x, base::DataVector& gradient) override {
    return fGradient(x, gradient);
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const override {
    clone =
        std::unique_ptr<ScalarFunctionGradient>(new WrapperScalarFunctionGradient(d, fGradient));
  }

 protected:
  /// function gradient to be wrapped
  FunctionGradientEvalType fGradient;
};
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUNCTION_SCALAR_WRAPPERSCALARFUNCTIONGRADIENT_HPP */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <functional>
#include <memory>
namespace sgpp {
namespace base {

/**
 * Implementation of ScalarFunction that wraps a std::function object.
 */
class WrapperScalarFunction : public ScalarFunction {
 public:
  typedef std::function<double(const DataVector&)> FunctionEvalType;

  /**
   * Constructor.
   *
   * @param d         dimension of the domain
   * @param f         function to be wrapped
   */
  WrapperScalarFunction(size_t d, FunctionEvalType f) : ScalarFunction(d), f(f) {}

  /**
   * Destructor.
   */
  ~WrapperScalarFunction() override {}

  /**
   * @param x     evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @return      \f$f(\vec{x})\f$
   */
  inline double eval(const DataVector& x) override { return f(x); }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<ScalarFunction>& clone) const override {
    clone = std::unique_ptr<ScalarFunction>(new WrapperScalarFunction(d, f));
  }

 protected:
  /// function to be wrapped
  FunctionEvalType f;
};
}  // namespace base
}  // namespace sgpp

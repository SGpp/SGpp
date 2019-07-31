// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <functional>
#include <memory>
namespace sgpp {
namespace base {

/**
 * Implementation of VectorFunction that wraps a std::function object.
 */
class WrapperVectorFunction : public VectorFunction {
 public:
  typedef std::function<void(const DataVector&, DataVector&)> FunctionEvalType;

  /**
   * Constructor.
   *
   * @param d         dimension of the domain
   * @param m         number of components
   * @param f         function to be wrapped
   */
  WrapperVectorFunction(size_t d, size_t m, FunctionEvalType f) : VectorFunction(d, m), f(f) {}

  /**
   * Destructor.
   */
  ~WrapperVectorFunction() override {}

  /**
   * @param[in]  x      evaluation point \f$\vec{x} \in [0, 1]^d\f$
   * @param[out] value  \f$g(\vec{x})\f$
   */
  inline void eval(const DataVector& x, DataVector& value) override { f(x, value); }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<VectorFunction>& clone) const override {
    clone = std::unique_ptr<VectorFunction>(new WrapperVectorFunction(d, m, f));
  }

 protected:
  /// function to be wrapped
  FunctionEvalType f;
};
}  // namespace base
}  // namespace sgpp

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/function/vector/WrapperVectorFunction.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Singleton containing an empty implementation of VectorFunction.
 * This is intended as a fill-in for ConstrainedOptimizer, if
 * only equality or inequality constraints are supported.
 */
class EmptyVectorFunction {
 public:
  static WrapperVectorFunction& getInstance();

 private:
  EmptyVectorFunction() {}
  EmptyVectorFunction(const EmptyVectorFunction&) = delete;
  void operator=(const EmptyVectorFunction&) = delete;
};
}  // namespace base
}  // namespace sgpp

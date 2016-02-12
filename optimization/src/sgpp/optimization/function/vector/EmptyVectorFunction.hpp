// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTION_HPP
#define SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTION_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/vector/WrapperVectorFunction.hpp>

namespace SGPP {
namespace optimization {

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
}  // namespace optimization
}  // namespace SGPP

#endif /* SGPP_OPTIMIZATION_FUNCTION_VECTOR_EMPTYVECTORFUNCTION_HPP */

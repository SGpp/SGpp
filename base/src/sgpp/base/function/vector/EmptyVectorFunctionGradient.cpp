// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/function/vector/EmptyVectorFunctionGradient.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

WrapperVectorFunctionGradient& EmptyVectorFunctionGradient::getInstance() {
  static WrapperVectorFunctionGradient wrapperVectorFunctionGradient(
      0, 0, [](const base::DataVector& x, base::DataVector& value, base::DataMatrix& gradient) {});
  return wrapperVectorFunctionGradient;
}
}  // namespace base
}  // namespace sgpp

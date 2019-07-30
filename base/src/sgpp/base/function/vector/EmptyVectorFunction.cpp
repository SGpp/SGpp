// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/function/vector/EmptyVectorFunction.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

WrapperVectorFunction& EmptyVectorFunction::getInstance() {
  static WrapperVectorFunction wrapperVectorFunction(
      0, 0, [](const base::DataVector& x, base::DataVector& value) {});
  return wrapperVectorFunction;
}
}  // namespace base
}  // namespace sgpp

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace combigrid {
/**
 * This is a helper class used internally to wrap python functions into GeneralFunction objects.
 * F should be a GeneralFunction<...>
 */
template <typename F>
class GeneralFunctionDirector {
 public:
  virtual ~GeneralFunctionDirector() {}

  virtual typename F::output_type eval(typename F::input_type vec) = 0;

  virtual F toFunction();
};

/**
 * This is a protected helper class used internally to wrap python functions into MultiFunction
 * objects.
 */
template <typename F>
class GFDirWrapper {
 public:
  GeneralFunctionDirector<F> *ptr;

  explicit GFDirWrapper(GeneralFunctionDirector<F> *ptr) : ptr(ptr) {}

  ~GFDirWrapper() { delete ptr; }
};

template <typename F>
F GeneralFunctionDirector<F>::toFunction() {
  auto wrapperPtr = std::make_shared<GFDirWrapper<F>>(this);
  return F([wrapperPtr](typename F::input_type in) ->
           typename F::output_type { return wrapperPtr->ptr->eval(in); });
}

} /* namespace combigrid */
} /* namespace sgpp */

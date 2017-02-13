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
template <typename T>
class GeneralFunctionDirector;

/**
 * This is a helper class used internally to wrap python functions into GeneralFunction objects.
 * F should be a GeneralFunction<...>
 */
template <typename Out, typename In>
class GeneralFunctionDirector<Out(In)> {
 public:
  virtual ~GeneralFunctionDirector() {}

  virtual Out eval(In vec) = 0;

  virtual GeneralFunction<Out(In)> toFunction();
};

template <typename Out>
class GeneralFunctionDirector<Out()> {
 public:
  virtual ~GeneralFunctionDirector() {}

  virtual Out eval() = 0;

  virtual GeneralFunction<Out()> toFunction();
};

template <typename T>
class GFDirWrapper;

/**
 * This is a protected helper class used internally to wrap python functions into MultiFunction
 * objects.
 */
template <typename Out, typename In>
class GFDirWrapper<Out(In)> {
 public:
  GeneralFunctionDirector<Out(In)> *ptr;

  explicit GFDirWrapper(GeneralFunctionDirector<Out(In)> *ptr) : ptr(ptr) {}

  ~GFDirWrapper() { delete ptr; }
};

template <typename Out>
class GFDirWrapper<Out()> {
 public:
  GeneralFunctionDirector<Out()> *ptr;

  explicit GFDirWrapper(GeneralFunctionDirector<Out()> *ptr) : ptr(ptr) {}

  ~GFDirWrapper() { delete ptr; }
};

template <typename Out, typename In>
GeneralFunction<Out(In)> GeneralFunctionDirector<Out(In)>::toFunction() {
  auto wrapperPtr = std::make_shared<GFDirWrapper<Out(In)>>(this);
  return F([wrapperPtr](In in) -> Out { return wrapperPtr->ptr->eval(in); });
}

template <typename Out>
GeneralFunction<Out()> GeneralFunctionDirector<Out()>::toFunction() {
  auto wrapperPtr = std::make_shared<GFDirWrapper<Out()>>(this);
  return F([wrapperPtr]() -> Out { return wrapperPtr->ptr->eval(); });
}

} /* namespace combigrid */
} /* namespace sgpp */

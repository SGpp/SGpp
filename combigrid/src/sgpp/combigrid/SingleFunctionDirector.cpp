// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/SingleFunctionDirector.hpp>

#include <memory>

namespace sgpp {
namespace combigrid {

/**
 * This is a protected helper class used internally to wrap python functions into MultiFunction
 * objects.
 */
class SFDirWrapper {
 public:
  SingleFunctionDirector *ptr;

  explicit SFDirWrapper(SingleFunctionDirector *ptr) : ptr(ptr) {}

  ~SFDirWrapper() { delete ptr; }
};

SingleFunctionDirector::~SingleFunctionDirector() {}

SingleFunction SingleFunctionDirector::toSingleFunction() {
  auto wrapperPtr = std::make_shared<SFDirWrapper>(this);
  return SingleFunction([wrapperPtr](double x) -> double { return wrapperPtr->ptr->eval(x); });
}

} /* namespace combigrid */
} /* namespace sgpp */

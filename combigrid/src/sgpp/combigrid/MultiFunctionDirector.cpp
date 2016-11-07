// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/MultiFunctionDirector.hpp>

#include <memory>

namespace sgpp {
namespace combigrid {

/**
 * This is a protected helper class used internally to wrap python functions into MultiFunction
 * objects.
 */
class MFDirWrapper {
 public:
  MultiFunctionDirector *ptr;

  explicit MFDirWrapper(MultiFunctionDirector *ptr) : ptr(ptr) {}

  ~MFDirWrapper() { delete ptr; }
};

MultiFunctionDirector::~MultiFunctionDirector() {}

MultiFunction MultiFunctionDirector::toMultiFunction() {
  auto wrapperPtr = std::make_shared<MFDirWrapper>(this);
  return MultiFunction(
      [wrapperPtr](base::DataVector const &vec) -> double { return wrapperPtr->ptr->eval(vec); });
}

} /* namespace combigrid */
} /* namespace sgpp */

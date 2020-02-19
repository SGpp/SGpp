// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <CL/cl.h>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>

namespace sgpp {
namespace base {

class OCLClonedBuffer {
 public:
  std::shared_ptr<OCLManager> manager;
  bool initialized;
  cl_mem* bufferList;
  size_t sizeofType;
  size_t elements;

 public:
  explicit OCLClonedBuffer(std::shared_ptr<OCLManager> manager);

  ~OCLClonedBuffer();

  bool isInitialized();

  cl_mem* getBuffer(size_t deviceNumber);

  void writeToBuffer(void* hostData, size_t* offsets = nullptr);

  void readFromBuffer(void* hostData, size_t* offsets = nullptr,
                      size_t* ranges = nullptr);

  void initializeBuffer(void* initialValues, size_t sizeofType,
                        size_t elements);

  void freeBuffer();
};

}  // namespace base
}  // namespace sgpp

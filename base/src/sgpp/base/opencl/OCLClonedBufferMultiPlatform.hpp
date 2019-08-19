// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <CL/cl.h>
#include <map>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>

namespace sgpp {
namespace base {

class OCLClonedBufferMultiPlatform {
 public:
  std::shared_ptr<OCLManagerMultiPlatform> manager;
  bool initialized;
  std::map<cl_platform_id, cl_mem*> platformBufferList;
  size_t sizeofType;
  size_t elements;

 public:
  explicit OCLClonedBufferMultiPlatform(
      std::shared_ptr<OCLManagerMultiPlatform> manager);

  ~OCLClonedBufferMultiPlatform();

  bool isInitialized();

  cl_mem* getBuffer(cl_platform_id platformId, size_t deviceIndex);

  void writeToBuffer(void* hostData, size_t* offsets = nullptr);

  void readFromBuffer(void* hostData, size_t* offsets = nullptr,
                      size_t* ranges = nullptr);

  void initializeBuffer(void* initialValues, size_t sizeofType,
                        size_t elements);

  void freeBuffer();
};

}  // namespace base
}  // namespace sgpp

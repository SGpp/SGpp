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

// copies the whole buffer on all devices, retrieves only the part
// that was worked on on a specific device

class OCLStretchedBufferMultiPlatform {
 private:
  std::shared_ptr<OCLManagerMultiPlatform> manager;
  bool initialized;
  std::map<cl_platform_id, cl_mem*> platformBufferList;
  size_t sizeofType;
  size_t elements;

  bool isMappedMemory;
  std::map<cl_platform_id, cl_mem> hostBuffer;
  std::map<cl_platform_id, void*> mappedHostBuffer;

 public:
  OCLStretchedBufferMultiPlatform(
      std::shared_ptr<OCLManagerMultiPlatform> manager);

  ~OCLStretchedBufferMultiPlatform();

  bool isInitialized();

  cl_mem* getBuffer(cl_platform_id platformId, size_t deviceNumber);

  // only for mapped buffer
  void* getMappedHostBuffer(cl_platform_id platformId);

  void copyToOtherHostBuffers(cl_platform_id originPlatformId);

  void readFromBuffer(std::map<cl_platform_id, size_t*> indexStart,
                      std::map<cl_platform_id, size_t*> indexEnd);

  void combineBuffer(std::map<cl_platform_id, size_t*> indexStart,
                     std::map<cl_platform_id, size_t*> indexEnd,
                     cl_platform_id platformId);

  void writeToBuffer();

  void initializeBuffer(size_t sizeofType, size_t elements);

  void freeBuffer();
};

}  // namespace base
}  // namespace sgpp

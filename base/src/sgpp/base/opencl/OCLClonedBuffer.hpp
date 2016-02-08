/*
 * OCLMemory.hpp
 *
 *  Created on: Mar 27, 2015
 *      Author: pfandedd
 */

#pragma once

#include <CL/cl.h>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/opencl/OCLManager.hpp>

namespace SGPP {
namespace base {

class OCLClonedBuffer {
 public:
  std::shared_ptr<OCLManager> manager;
  bool initialized;
  cl_mem* bufferList;
  size_t sizeofType;
  size_t elements;

 public:

  OCLClonedBuffer(std::shared_ptr<OCLManager> manager);

  ~OCLClonedBuffer();

  bool isInitialized();

  cl_mem* getBuffer(size_t deviceNumber);

  void writeToBuffer(void* hostData, size_t* offsets = nullptr);

  void readFromBuffer(void* hostData, size_t* offsets = nullptr,
                      size_t* ranges = nullptr);

  void initializeBuffer(void* initialValues, size_t sizeofType, size_t elements);

  void freeBuffer();

};

}  // namespace base
}  // namespace SGPP



/*
 * OCLStretchedBuffer.hpp
 *
 *  Created on: Mar 31, 2015
 *      Author: pfandedd
 */

#pragma once

#include <CL/cl.h>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/opencl/OCLManager.hpp>

namespace SGPP {
namespace base {

//copies the whole buffer on all devices, retrieves only the part that was worked on on a specific device

class OCLStretchedBuffer {
 private:
  std::shared_ptr<OCLManager> manager;
  bool initialized;
  cl_mem* bufferList;
  size_t sizeofType;
  size_t elements;

  bool isMappedMemory;
  cl_mem hostBuffer;
  void* mappedHostBuffer;
 public:

  OCLStretchedBuffer(std::shared_ptr<OCLManager> manager);

  ~OCLStretchedBuffer();

  bool isInitialized();

  cl_mem* getBuffer(size_t deviceNumber);

  //only for mapped buffer
  void* getMappedHostBuffer();

  void readFromBuffer(size_t* indexStart, size_t* indexEnd);

  void writeToBuffer();

  void initializeBuffer(size_t sizeofType, size_t elements);

  void freeBuffer();
};

}  // namespace base
}  // namespace SGPP



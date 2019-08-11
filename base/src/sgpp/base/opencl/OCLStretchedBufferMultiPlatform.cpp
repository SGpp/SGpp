// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <cstring>  // required for "memmove"
#include <string>
#include <sstream>
#include <map>
#include <vector>

#include <sgpp/base/opencl/OCLStretchedBufferMultiPlatform.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

namespace sgpp {
namespace base {

OCLStretchedBufferMultiPlatform::OCLStretchedBufferMultiPlatform(
    std::shared_ptr<OCLManagerMultiPlatform> manager)
    : manager(manager) {
  initialized = false;
  sizeofType = 0;
  elements = 0;
  isMappedMemory = false;
}

OCLStretchedBufferMultiPlatform::~OCLStretchedBufferMultiPlatform() { this->freeBuffer(); }

bool OCLStretchedBufferMultiPlatform::isInitialized() { return this->initialized; }

cl_mem* OCLStretchedBufferMultiPlatform::getBuffer(cl_platform_id platformId, size_t deviceNumber) {
  return &(this->platformBufferList[platformId][deviceNumber]);
}

void OCLStretchedBufferMultiPlatform::initializeBuffer(size_t sizeofType, size_t elements) {
  cl_int err;

  // cl buffer that is allocated on the host, but cannot be directly accessed
  for (OCLPlatformWrapper& platform : manager->platforms) {
    cl_mem hostBuffer = clCreateBuffer(platform.context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                       sizeofType * elements, nullptr, &err);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Could not allocate host buffer! "
                     "Error code: " << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    this->hostBuffer[platform.platformId] = hostBuffer;

    cl_mem* bufferList = new cl_mem[platform.getDeviceCount()];

    for (size_t i = 0; i < platform.getDeviceCount(); i++) {
      bufferList[i] =
          clCreateBuffer(platform.context, CL_MEM_READ_WRITE, sizeofType * elements, nullptr, &err);

      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Could not allocate buffer! "
                       "Error code: " << err << std::endl;
        throw sgpp::base::operation_exception(errorString.str());
      }
    }

    this->platformBufferList[platform.platformId] = bufferList;

    void* hostPinnedMemory = clEnqueueMapBuffer(platform.commandQueues[0], hostBuffer, CL_TRUE,
                                                CL_MAP_READ | CL_MAP_WRITE, 0,
                                                sizeofType * elements, 0, nullptr, nullptr, &err);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Could not map pinned memory to host pointer! "
                     "Error code: " << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    this->mappedHostBuffer[platform.platformId] = hostPinnedMemory;
  }

  this->sizeofType = sizeofType;
  this->elements = elements;
  this->isMappedMemory = true;
  this->initialized = true;
}

void OCLStretchedBufferMultiPlatform::freeBuffer() {
  if (!this->initialized) {
    // permitted to allow for easy kernel resets
    return;
  }

  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    if (this->mappedHostBuffer[platform.platformId] == nullptr) {
      std::stringstream errorString;
      errorString << "OCL Error: OCLStretchedBufferMultiPlatform in "
                     "partially initialized state: mappedHostBuffer is null" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    if (this->hostBuffer[platform.platformId] == nullptr) {
      std::stringstream errorString;
      errorString << "OCL Error: OCLStretchedBufferMultiPlatform in "
                     "partially initialized state: hostBuffer is null" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    cl_int err =
        clEnqueueUnmapMemObject(platform.commandQueues[0], hostBuffer[platform.platformId],
                                this->mappedHostBuffer[platform.platformId], 0, nullptr, nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: OCLStretchedBuffer unmapping "
                     "memory not successful" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    if (platformBufferList[platform.platformId] == nullptr) {
      std::stringstream errorString;
      errorString << "OCL Error: OCLStretchedBufferMultiPlatform in "
                     "partially initialized state: platform buffer list is null" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    for (size_t i = 0; i < platform.getDeviceCount(); i++) {
      if (this->platformBufferList[platform.platformId][i] != nullptr) {
        clReleaseMemObject(this->platformBufferList[platform.platformId][i]);
        this->platformBufferList[platform.platformId][i] = nullptr;
      } else {
        std::stringstream errorString;
        errorString << "OCL Error: OCLStretchedBufferMultiPlatform in "
                       "partially initialized state: device buffer is null" << std::endl;
        throw sgpp::base::operation_exception(errorString.str());
      }
    }

    delete[] this->platformBufferList[platform.platformId];
    this->platformBufferList[platform.platformId] = nullptr;

    if (this->isMappedMemory) {
      if (this->hostBuffer[platform.platformId] != nullptr) {
        clReleaseMemObject(this->hostBuffer[platform.platformId]);
        this->hostBuffer[platform.platformId] = nullptr;
      } else {
        std::stringstream errorString;
        errorString << "OCL Error: OCLStretchedBufferMultiPlatform in "
                       "partially initialized state: host buffer is null" << std::endl;
        throw sgpp::base::operation_exception(errorString.str());
      }
    }
  }

  this->initialized = false;
}

void* OCLStretchedBufferMultiPlatform::getMappedHostBuffer(cl_platform_id platformId) {
  return this->mappedHostBuffer[platformId];
}

void OCLStretchedBufferMultiPlatform::copyToOtherHostBuffers(cl_platform_id originPlatformId) {
  if (!this->initialized) {
    std::stringstream errorString;
    errorString << "OCL Error: tried to \"copyToOtherHostBuffers\" "
                   "with uninitialized buffer" << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    if (platform.platformId == originPlatformId) {
      continue;
    }

    memmove(this->mappedHostBuffer[platform.platformId], this->mappedHostBuffer[originPlatformId],
            this->elements * this->sizeofType);
  }
}

void OCLStretchedBufferMultiPlatform::readFromBuffer(std::map<cl_platform_id, size_t*> indexStart,
                                                     std::map<cl_platform_id, size_t*> indexEnd) {
  cl_int err = CL_SUCCESS;
  //    cl_event* actionDone = new cl_event[this->manager->overallDeviceCount];
  std::map<cl_platform_id, std::vector<cl_event>> platformActionEvents;
  std::map<cl_platform_id, size_t> platformTransferringDevice;

  // read data back
  //    size_t devicesTransferring = 0;
  //    size_t actionIndex = 0;

  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    platformActionEvents[platform.platformId] = std::vector<cl_event>(platform.getDeviceCount());
    size_t devicesTransferring = 0;

    for (size_t i = 0; i < platform.getDeviceCount(); i++) {
      size_t indexEndDevice = indexEnd[platform.platformId][i];
      size_t indexStartDevice = indexStart[platform.platformId][i];
      size_t range = indexEndDevice - indexStartDevice;

      if (range != 0) {
        err = clEnqueueReadBuffer(platform.commandQueues[i],
                                  platformBufferList[platform.platformId][i], CL_FALSE,
                                  this->sizeofType * indexStartDevice, this->sizeofType * range,
                                  static_cast<char*>(this->mappedHostBuffer[platform.platformId]) +
                                      (this->sizeofType * indexStartDevice),
                                  0, nullptr, &(platformActionEvents[platform.platformId][i]));
        devicesTransferring += 1;
      } else {
        break;
      }

      if (err != CL_SUCCESS && range != 0) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to enqueue read buffer command! "
                       "Error code: " << err << std::endl;
        throw sgpp::base::operation_exception(errorString.str());
      }
    }

    platformTransferringDevice[platform.platformId] = devicesTransferring;
  }

  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    std::vector<cl_event> events = platformActionEvents[platform.platformId];
    clWaitForEvents((cl_uint)platformTransferringDevice[platform.platformId], events.data());
  }

  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    size_t devicesTransferring = platformTransferringDevice[platform.platformId];

    for (size_t i = 0; i < devicesTransferring; i++) {
      cl_event event = platformActionEvents[platform.platformId][i];
      clReleaseEvent(event);
    }
  }
}

void OCLStretchedBufferMultiPlatform::combineBuffer(std::map<cl_platform_id, size_t*> indexStart,
                                                    std::map<cl_platform_id, size_t*> indexEnd,
                                                    cl_platform_id platformId) {
  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    if (platform.platformId == platformId) {
      continue;
    }

    for (size_t deviceIndex = 0; deviceIndex < platform.getDeviceCount(); deviceIndex++) {
      size_t indexStartDevice = indexStart[platform.platformId][deviceIndex];
      size_t indexEndDevice = indexEnd[platform.platformId][deviceIndex];
      size_t range = indexEndDevice - indexStartDevice;
      memmove(static_cast<char*>(this->mappedHostBuffer[platformId]) +
                  indexStartDevice * this->sizeofType,
              static_cast<char*>(this->mappedHostBuffer[platform.platformId]) +
                  indexStartDevice * this->sizeofType,
              range * this->sizeofType);
    }
  }
}

void OCLStretchedBufferMultiPlatform::writeToBuffer() {
  cl_int err;

  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    for (size_t i = 0; i < platform.getDeviceCount(); i++) {
      err = clEnqueueWriteBuffer(platform.commandQueues[i],
                                 platformBufferList[platform.platformId][i], CL_TRUE, 0,
                                 this->sizeofType * this->elements,
                                 this->mappedHostBuffer[platform.platformId], 0, nullptr, nullptr);

      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to enqueue mapped write command! "
                       "Error code: " << err << std::endl;
        throw sgpp::base::operation_exception(errorString.str());
      }
    }
  }
}

}  // namespace base
}  // namespace sgpp

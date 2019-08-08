// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sstream>
#include <map>
#include <vector>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLClonedBufferMultiPlatform.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

namespace sgpp {
namespace base {

OCLClonedBufferMultiPlatform::OCLClonedBufferMultiPlatform(
    std::shared_ptr<OCLManagerMultiPlatform> manager)
    : manager(manager) {
  initialized = false;
  sizeofType = 0;
  elements = 0;
}

OCLClonedBufferMultiPlatform::~OCLClonedBufferMultiPlatform() { this->freeBuffer(); }

bool OCLClonedBufferMultiPlatform::isInitialized() { return this->initialized; }

cl_mem* OCLClonedBufferMultiPlatform::getBuffer(cl_platform_id platformId, size_t deviceIndex) {
  return &(this->platformBufferList[platformId][deviceIndex]);
}

void OCLClonedBufferMultiPlatform::writeToBuffer(void* hostData, size_t* offsets) {
  cl_int err;

  std::map<cl_platform_id, std::vector<cl_event>> platformActionEvents;
  // cl_event* actionDone = new cl_event[this->manager->overallDeviceCount];

  // size_t actionIndex = 0;
  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    platformActionEvents[platform.platformId] = std::vector<cl_event>(platform.getDeviceCount());

    for (size_t i = 0; i < platform.getDeviceCount(); i++) {
      if (offsets == nullptr) {
        err = clEnqueueWriteBuffer(platform.commandQueues[i],
                                   this->platformBufferList[platform.platformId][i], CL_FALSE, 0,
                                   this->sizeofType * this->elements, hostData, 0, nullptr,
                                   &(platformActionEvents[platform.platformId][i]));
      } else {
        err = clEnqueueWriteBuffer(
            platform.commandQueues[i], this->platformBufferList[platform.platformId][i], CL_FALSE,
            this->sizeofType * offsets[i], this->sizeofType * this->elements, hostData, 0, nullptr,
            &(platformActionEvents[platform.platformId][i]));
      }

      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to enqueue write buffer command! "
                       "Error code: " << err << std::endl;
        throw sgpp::base::operation_exception(errorString.str());
      }

      // actionIndex += 1;
    }
  }

  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    std::vector<cl_event>& events = platformActionEvents[platform.platformId];
    clWaitForEvents(static_cast<cl_uint>(platform.getDeviceCount()), events.data());
  }

  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    for (size_t i = 0; i < platform.getDeviceCount(); i++) {
      cl_event event = platformActionEvents[platform.platformId][i];
      clReleaseEvent(event);
    }
  }
  platformActionEvents.clear();
}

void OCLClonedBufferMultiPlatform::readFromBuffer(void* hostData, size_t* offsets, size_t* ranges) {
  cl_int err;

  //    cl_event* actionDone = new cl_event[this->manager->overallDeviceCount];
  std::map<cl_platform_id, std::vector<cl_event>> platformActionEvents;

  size_t actionIndex = 0;

  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    platformActionEvents[platform.platformId] = std::vector<cl_event>(platform.getDeviceCount());

    // read data back
    for (size_t i = 0; i < platform.getDeviceCount(); i++) {
      if (offsets == nullptr) {
        err = clEnqueueReadBuffer(platform.commandQueues[i],
                                  this->platformBufferList[platform.platformId][i], CL_FALSE, 0,
                                  this->sizeofType * this->elements, hostData, 0, nullptr,
                                  &(platformActionEvents[platform.platformId][i]));
      } else {
        err = clEnqueueReadBuffer(platform.commandQueues[i],
                                  this->platformBufferList[platform.platformId][i], CL_FALSE,
                                  this->sizeofType * offsets[i], this->sizeofType * ranges[i],
                                  static_cast<char*>(hostData) + (this->sizeofType * offsets[i]), 0,
                                  nullptr, &(platformActionEvents[platform.platformId][i]));
      }

      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to enqueue read buffer command! "
                       "Error code: " << err << std::endl;
        throw sgpp::base::operation_exception(errorString.str());
      }

      actionIndex += 1;
    }
  }

  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    std::vector<cl_event>& events = platformActionEvents[platform.platformId];
    clWaitForEvents(static_cast<cl_uint>(platform.getDeviceCount()), events.data());
  }

  for (OCLPlatformWrapper& platform : this->manager->platforms) {
    for (size_t i = 0; i < platform.getDeviceCount(); i++) {
      cl_event event = platformActionEvents[platform.platformId][i];
      clReleaseEvent(event);
    }
  }
  platformActionEvents.clear();
}

void OCLClonedBufferMultiPlatform::initializeBuffer(void* initialValues, size_t sizeofType,
                                                    size_t elements) {
  cl_int err;

  for (OCLPlatformWrapper& platform : manager->platforms) {
    cl_mem* bufferList = new cl_mem[platform.getDeviceCount()];

    for (size_t i = 0; i < platform.getDeviceCount(); i++) {
      if (initialValues != nullptr) {
        bufferList[i] = clCreateBuffer(platform.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                       sizeofType * elements, initialValues, &err);
      } else {
        bufferList[i] = clCreateBuffer(platform.context, CL_MEM_READ_ONLY, sizeofType * elements,
                                       nullptr, &err);
      }

      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Could not allocate buffer! "
                       "Error code: " << err << std::endl;
        throw sgpp::base::operation_exception(errorString.str());
      }
    }

    this->platformBufferList[platform.platformId] = bufferList;
  }

  this->sizeofType = sizeofType;
  this->elements = elements;
  this->initialized = true;
}

void OCLClonedBufferMultiPlatform::freeBuffer() {
  if (!this->initialized) {
    // permitted to allow for easy kernel resets
    return;
  }

  for (OCLPlatformWrapper& platform : manager->platforms) {
    if (this->platformBufferList[platform.platformId] == nullptr) {
      std::stringstream errorString;
      errorString << "OCL Error: OCLClonedBufferMultiPlatform in partially "
                     "initialized state: platform buffer list is null" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    cl_mem* bufferList = this->platformBufferList[platform.platformId];

    for (size_t i = 0; i < platform.getDeviceCount(); i++) {
      if (bufferList[i] != nullptr) {
        clReleaseMemObject(bufferList[i]);
        bufferList[i] = nullptr;
      } else {
        std::stringstream errorString;
        errorString << "OCL Error: OCLClonedBufferMultiPlatform in partially "
                       "initialized state: device buffer is null" << std::endl;
        throw sgpp::base::operation_exception(errorString.str());
      }
    }

    delete[] this->platformBufferList[platform.platformId];
    this->platformBufferList[platform.platformId] = nullptr;
  }

  this->initialized = false;
}

}  // namespace base
}  // namespace sgpp

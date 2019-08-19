// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <CL/cl.h>
#include <vector>

#include <sgpp/base/opencl/OCLDevice.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

template <typename T> class OCLBufferWrapperSD {
private:
  std::shared_ptr<OCLDevice> device;
  bool initialized;
  cl_mem buffer;
  //    size_t sizeofType;
  size_t elements;
  std::vector<T> hostData;

public:
  explicit OCLBufferWrapperSD(std::shared_ptr<base::OCLDevice> device)
      : device(device), initialized(false), buffer(nullptr), elements(0) {}

  ~OCLBufferWrapperSD() { this->freeBuffer(); }

  bool isInitialized() { return this->initialized; }

  cl_mem *getBuffer() {
    if (!this->initialized) {
      std::stringstream errorString;
      errorString << "OCL Error: Buffer not initialized: " << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
    return &this->buffer;
  }

  size_t size() {
    if (!this->initialized) {
      std::stringstream errorString;
      errorString << "OCL Error: Buffer not initialized: " << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
    return this->elements;
  }

  std::vector<T> &getHostPointer() {
    if (!this->initialized) {
      std::stringstream errorString;
      errorString << "OCL Error: Buffer not initialized: " << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
    return hostData;
  }

  void writeToBuffer() {
    if (!this->initialized) {
      std::stringstream errorString;
      errorString << "OCL Error: Buffer not initialized: " << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    cl_int err;

    err = clEnqueueWriteBuffer(device->commandQueue, this->buffer, CL_TRUE, 0,
                               sizeof(T) * this->elements, hostData.data(), 0,
                               nullptr, nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString
          << "OCL Error: Failed to enqueue write buffer command! Error code: "
          << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  void readFromBuffer() {
    if (!this->initialized) {
      std::stringstream errorString;
      errorString << "OCL Error: Buffer not initialized: " << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    cl_int err;

    err = clEnqueueReadBuffer(device->commandQueue, this->buffer, CL_TRUE, 0,
                              sizeof(T) * this->elements,
                              static_cast<void *>(hostData.data()), 0, nullptr,
                              nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString
          << "OCL Error: Failed to enqueue read buffer command! Error code: "
          << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  void initializeBuffer(size_t elements) {
    cl_int err;

    hostData.resize(elements);

    //        if (initialValues.size() != 0) {
    //            if (initialValues.size() != elements) {
    //                std::stringstream errorString;
    //                errorString << "OCL Error: Size of initial values vector
    //                does not match specified size! Error code: "
    //                        << err << std::endl;
    //                throw sgpp::base::operation_exception(errorString.str());
    //            }
    //            this->buffer = clCreateBuffer(device->context,
    //            CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(T) *
    //            elements, initialValues.data(), &err);
    //        } else {
    this->buffer = clCreateBuffer(device->context, CL_MEM_READ_WRITE,
                                  sizeof(T) * elements, nullptr, &err);
    //        }

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Could not allocate buffer! Error code: " << err
                  << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    this->elements = elements;
    this->initialized = true;
  }

  void freeBuffer() {
    if (!this->initialized) {
      // permitted to allow for easy kernel resets
      return;
    }

    if (buffer != nullptr) {
      clReleaseMemObject(buffer);
      buffer = nullptr;
    } else {
      std::stringstream errorString;
      errorString << "OCL Error: could not free OCLClonedBufferSD" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
    this->initialized = false;
    this->elements = 0;
  }

  void intializeTo(std::vector<T> &hostBuffer, size_t dim, size_t offsetStart,
                   size_t offsetEnd, bool storeStructOfArrays = false) {
    size_t range = offsetEnd - offsetStart;
    size_t totalElements = range * dim;

    if (!this->isInitialized()) {
      this->initializeBuffer(totalElements);
    } else if (totalElements != this->size()) {
      this->freeBuffer();
      this->initializeBuffer(totalElements);
    }

    size_t dataPoints = hostBuffer.size() / dim;

    if (dataPoints < range) {
      std::stringstream errorString;
      errorString << "OCL Error: initializeTo requires hostbuffer size >= "
                     "device buffer size"
                  << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    std::vector<T> &deviceDataHost = this->getHostPointer();

    // memory layout: AAABBBCCC -> struct of arrays (SOA)
    if (storeStructOfArrays) {
      for (size_t d = 0; d < dim; d++) {
        size_t deviceDataIndex = 0;
        for (size_t i = offsetStart; i < offsetEnd; i++) {
          deviceDataHost[d * range + deviceDataIndex] =
              hostBuffer[d * dataPoints + i];
          deviceDataIndex += 1;
        }
      }
    } else {
      // memory layout: ABCABCABC -> AOS
      size_t deviceDataIndex = 0;
      for (size_t i = offsetStart; i < offsetEnd; i++) {
        for (size_t d = 0; d < dim; d++) {
          deviceDataHost[deviceDataIndex * dim + d] = hostBuffer[i * dim + d];
        }
        deviceDataIndex += 1;
      }
    }

    this->writeToBuffer();
  }
};
} // namespace base
} // namespace sgpp

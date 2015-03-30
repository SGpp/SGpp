/*
 * OCLMemory.cpp
 *
 *  Created on: Mar 30, 2015
 *      Author: pfandedd
 */

#include <sstream>

#include <sgpp/globaldef.hpp>

#include "OCLMemory.hpp"

#include <sgpp/base/exception/operation_exception.hpp>

namespace SGPP {
namespace datadriven {

OCLMemory::OCLMemory(OCLManager &manager) :
    manager(manager) {
  initialized = false;
  bufferList = nullptr;
  sizeofType = 0;
  elements = 0;
  isMappedMemory = false;
  hostBuffer = 0;
  mappedHostBuffer = nullptr;
}

bool OCLMemory::isInitialized() {
  return this->initialized;
}

cl_mem *OCLMemory::getBuffer(size_t deviceNumber) {
  return &(this->bufferList[deviceNumber]);
}

//TODO: current multidevice strategy: allocate everything everywere, use only range specified for device

void OCLMemory::writeToBuffer(void *hostData, size_t *offsets) {
  cl_int err;
  cl_event* actionDone = new cl_event[this->manager.num_devices];

  for (size_t i = 0; i < this->manager.num_devices; i++) {
    if (offsets == nullptr) {
      err = clEnqueueWriteBuffer(this->manager.command_queue[i], this->bufferList[i],
      CL_FALSE, 0, this->sizeofType * this->elements, hostData, 0, nullptr, nullptr);
    } else {
      err = clEnqueueWriteBuffer(this->manager.command_queue[i], this->bufferList[i],
      CL_FALSE, this->sizeofType * offsets[i], this->sizeofType * this->elements, hostData, 0, nullptr, nullptr);
    }

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to enqueue write buffer command! Error code: " << err << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }
  }

  clWaitForEvents((cl_uint) this->manager.num_devices, actionDone);
//  for (size_t i = 0; i < this->manager.num_devices; i++) {
//    if (indexEnd[i] > indexStart[i]) {
//        clReleaseEvent(actionDone[i]);
//    }
//  }
}

void OCLMemory::readFromBuffer(void *hostData, size_t *offsets, size_t *ranges) {
  cl_int err;
  cl_event* actionDone = new cl_event[this->manager.num_devices];

  // read data back
  for (size_t i = 0; i < this->manager.num_devices; i++) {
    if (offsets == nullptr) {
      err = clEnqueueReadBuffer(this->manager.command_queue[i], this->bufferList[i],
      CL_FALSE, 0, this->sizeofType * this->elements, hostData, 0, nullptr, &(actionDone[i]));
    } else {
      err = clEnqueueReadBuffer(this->manager.command_queue[i], this->bufferList[i],
      CL_FALSE, this->sizeofType * offsets[i], this->sizeofType * ranges[i],
          static_cast<char *>(hostData) + (this->sizeofType * offsets[i]), 0, nullptr, &(actionDone[i]));
    }

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to enqueue read buffer command! Error code: " << err << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }
  }
  clWaitForEvents((cl_uint) this->manager.num_devices, actionDone);
}

//TODO: read/write-flags

void OCLMemory::initializeBuffer(void *initialValues, size_t sizeofType, size_t elements) {
  cl_int err;
  cl_mem *bufferList = new cl_mem[manager.num_devices];
  for (size_t i = 0; i < manager.num_devices; i++) {
    if (initialValues != nullptr) {
      bufferList[i] = clCreateBuffer(manager.context,
      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeofType * elements, initialValues, &err);
    } else {
      bufferList[i] = clCreateBuffer(manager.context,
      CL_MEM_READ_ONLY, sizeofType * elements, nullptr, &err);
    }
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Could not allocate buffer! Error code: " << err << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }
  }
  this->bufferList = bufferList;
  this->sizeofType = sizeofType;
  this->elements = elements;
  this->isMappedMemory = false;
  this->hostBuffer = nullptr;
  this->mappedHostBuffer = nullptr;
}

void OCLMemory::initializeMappedBuffer(size_t sizeofType, size_t elements) {
  cl_int err;
  //cl buffer that is allocated on the host, but cannot be directly accessed
  cl_mem hostBuffer = clCreateBuffer(manager.context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeofType * elements,
      nullptr, &err);
  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Could not allocate host buffer! Error code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  cl_mem *bufferList = new cl_mem[manager.num_devices];
  for (size_t i = 0; i < manager.num_devices; i++) {
    bufferList[i] = clCreateBuffer(manager.context, CL_MEM_READ_WRITE, sizeofType * elements, nullptr, &err);
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Could not allocate buffer! Error code: " << err << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }
  }

  //TODO: why command queue 0?
  void *hostPinnedMemory = clEnqueueMapBuffer(manager.command_queue[0], hostBuffer, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE,
      0, sizeofType * elements, 0, nullptr, nullptr, &err);
  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Could not map pinned memory to host pointer! Error code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  this->bufferList = bufferList;
  this->sizeofType = sizeofType;
  this->elements = elements;
  this->isMappedMemory = true;
  this->hostBuffer = hostBuffer;
  this->mappedHostBuffer = hostPinnedMemory;
}

void OCLMemory::freeBuffer() {
  if (!this->initialized) {
    return;
  }

  for (size_t i = 0; i < this->manager.num_devices; i++) {
    if (this->bufferList[i] != nullptr) {
      clReleaseMemObject(this->bufferList[i]);
      this->bufferList[i] = nullptr;
    }
  }
  delete[] this->bufferList;
  if (this->isMappedMemory) {
    clReleaseMemObject(this->hostBuffer);
  }
  this->initialized = false;
}

void *OCLMemory::getMappedHostBuffer() {
  return this->mappedHostBuffer;
}

void OCLMemory::readMappedBuffer(size_t *indexStart, size_t *indexEnd) {
  if (indexStart == nullptr || indexEnd == nullptr) {
    std::stringstream errorString;
    errorString << "OCL Error: reading mapped buffer failed, inconsistent arguments used (see documentation)"
        << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }
  cl_int err = CL_SUCCESS;
  cl_event* actionDone = new cl_event[this->manager.num_devices];

  // read data back
  size_t devicesTransferring = 0;
  for (size_t i = 0; i < this->manager.num_devices; i++) {
    size_t range = indexEnd[i] - indexStart[i];
    if (range != 0) {
      err = clEnqueueReadBuffer(this->manager.command_queue[i], this->bufferList[i],
      CL_FALSE, this->sizeofType * indexStart[i], this->sizeofType * range,
          static_cast<char *>(this->mappedHostBuffer) + (this->sizeofType * indexStart[i]), 0, nullptr,
          &(actionDone[i]));
      devicesTransferring += 1;
    }

    if (err != CL_SUCCESS && range != 0) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to enqueue read buffer command! Error code: " << err << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }
  }

  if (devicesTransferring > 0) {
    clWaitForEvents((cl_uint) this->manager.num_devices, actionDone);
  }

  for (size_t i = 0; i < this->manager.num_devices; i++) {
    if (indexEnd[i] > indexStart[i]) {
        clReleaseEvent(actionDone[i]);
    }
  }
}

void OCLMemory::writeMappedBuffer() {
  cl_int err;
  for (size_t i = 0; i < manager.num_devices; i++) {
    err = clEnqueueWriteBuffer(manager.command_queue[i], bufferList[i],
    CL_TRUE, 0, this->sizeofType * this->elements, this->mappedHostBuffer, 0, nullptr, nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to enqueue mapped write command! Error code: " << err << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }
  }
}

}
}


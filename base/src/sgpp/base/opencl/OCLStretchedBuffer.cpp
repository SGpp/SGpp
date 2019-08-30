// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sstream>

#include <sgpp/base/opencl/OCLStretchedBuffer.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

namespace sgpp {
namespace base {

OCLStretchedBuffer::OCLStretchedBuffer(std::shared_ptr<OCLManager> manager) : manager(manager) {
  initialized = false;
  bufferList = nullptr;
  sizeofType = 0;
  elements = 0;
  isMappedMemory = false;
  hostBuffer = 0;
  mappedHostBuffer = nullptr;
}

OCLStretchedBuffer::~OCLStretchedBuffer() { this->freeBuffer(); }

bool OCLStretchedBuffer::isInitialized() { return this->initialized; }

cl_mem* OCLStretchedBuffer::getBuffer(size_t deviceNumber) {
  return &(this->bufferList[deviceNumber]);
}

void OCLStretchedBuffer::initializeBuffer(size_t sizeofType, size_t elements) {
  cl_int err;
  // cl buffer that is allocated on the host, but cannot be directly accessed
  cl_mem hostBuffer = clCreateBuffer(manager->context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                     sizeofType * elements, nullptr, &err);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Could not allocate host buffer! "
                   "Error code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  cl_mem* bufferList = new cl_mem[manager->num_devices];

  for (size_t i = 0; i < manager->num_devices; i++) {
    bufferList[i] =
        clCreateBuffer(manager->context, CL_MEM_READ_WRITE, sizeofType * elements, nullptr, &err);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Could not allocate buffer! "
                     "Error code: " << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  void* hostPinnedMemory =
      clEnqueueMapBuffer(manager->command_queue[0], hostBuffer, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE,
                         0, sizeofType * elements, 0, nullptr, nullptr, &err);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Could not map pinned memory to host pointer! "
                   "Error code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  this->bufferList = bufferList;
  this->sizeofType = sizeofType;
  this->elements = elements;
  this->isMappedMemory = true;
  this->hostBuffer = hostBuffer;
  this->mappedHostBuffer = hostPinnedMemory;
  this->initialized = true;
}

void OCLStretchedBuffer::freeBuffer() {
  if (!this->initialized) {
    return;
  }

  if (this->mappedHostBuffer == nullptr) {
    std::stringstream errorString;
    errorString << "OCL Error: OCLStretchedBuffer in partially initialized "
                   "state: mappedHostBuffer is null" << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  if (this->hostBuffer == nullptr) {
    std::stringstream errorString;
    errorString << "OCL Error: OCLStretchedBuffer in partially initialized "
                   "state: hostBuffer is null" << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  cl_int err = clEnqueueUnmapMemObject(manager->command_queue[0], hostBuffer,
                                       this->mappedHostBuffer, 0, nullptr, nullptr);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: OCLStretchedBuffer unmapping memory "
                   "not successful" << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  if (this->bufferList == nullptr) {
    std::stringstream errorString;
    errorString << "OCL Error: OCLStretchedBuffer in partially "
                   "initialized state: buffer list is null" << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  for (size_t i = 0; i < this->manager->num_devices; i++) {
    if (this->bufferList[i] != nullptr) {
      clReleaseMemObject(this->bufferList[i]);
      this->bufferList[i] = nullptr;
    } else {
      std::stringstream errorString;
      errorString << "OCL Error: OCLStretchedBuffer in partially "
                     "initialized state: device buffer is null" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  delete[] this->bufferList;
  this->bufferList = nullptr;

  if (this->isMappedMemory) {
    if (this->hostBuffer != nullptr) {
      clReleaseMemObject(this->hostBuffer);
      this->hostBuffer = nullptr;
    } else {
      std::stringstream errorString;
      errorString << "OCL Error: OCLStretchedBuffer in partially "
                     "initialized state: host buffer is null" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  this->initialized = false;
}

void* OCLStretchedBuffer::getMappedHostBuffer() { return this->mappedHostBuffer; }

void OCLStretchedBuffer::readFromBuffer(size_t* indexStart, size_t* indexEnd) {
  if (indexStart == nullptr || indexEnd == nullptr) {
    std::stringstream errorString;
    errorString << "OCL Error: reading mapped buffer failed, "
                   "inconsistent arguments used (see documentation)" << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  cl_int err = CL_SUCCESS;
  cl_event* actionDone = new cl_event[this->manager->num_devices];

  // read data back
  size_t devicesTransferring = 0;

  for (size_t i = 0; i < this->manager->num_devices; i++) {
    size_t range = indexEnd[i] - indexStart[i];

    if (range != 0) {
      err = clEnqueueReadBuffer(
          this->manager->command_queue[i], this->bufferList[i], CL_FALSE,
          this->sizeofType * indexStart[i], this->sizeofType * range,
          static_cast<char*>(this->mappedHostBuffer) + (this->sizeofType * indexStart[i]), 0,
          nullptr, &(actionDone[i]));
      devicesTransferring += 1;
    }

    if (err != CL_SUCCESS && range != 0) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to enqueue read buffer command! "
                     "Error code: " << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  if (devicesTransferring > 0) {
    clWaitForEvents((cl_uint)devicesTransferring, actionDone);
  }

  for (size_t i = 0; i < this->manager->num_devices; i++) {
    if (indexEnd[i] > indexStart[i]) {
      clReleaseEvent(actionDone[i]);
    }
  }
}

void OCLStretchedBuffer::writeToBuffer() {
  cl_int err;

  for (size_t i = 0; i < manager->num_devices; i++) {
    err = clEnqueueWriteBuffer(manager->command_queue[i], bufferList[i], CL_TRUE, 0,
                               this->sizeofType * this->elements, this->mappedHostBuffer, 0,
                               nullptr, nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to enqueue mapped write command! "
                     "Error code: " << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }
}

}  // namespace base
}  // namespace sgpp

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sstream>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLClonedBuffer.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

namespace sgpp {
namespace base {

OCLClonedBuffer::OCLClonedBuffer(std::shared_ptr<OCLManager> manager) : manager(manager) {
  initialized = false;
  bufferList = nullptr;
  sizeofType = 0;
  elements = 0;
}

OCLClonedBuffer::~OCLClonedBuffer() { this->freeBuffer(); }

bool OCLClonedBuffer::isInitialized() { return this->initialized; }

cl_mem* OCLClonedBuffer::getBuffer(size_t deviceNumber) {
  return &(this->bufferList[deviceNumber]);
}

void OCLClonedBuffer::writeToBuffer(void* hostData, size_t* offsets) {
  cl_int err;
  cl_event* actionDone = new cl_event[this->manager->num_devices];

  for (size_t i = 0; i < this->manager->num_devices; i++) {
    if (offsets == nullptr) {
      err = clEnqueueWriteBuffer(this->manager->command_queue[i], this->bufferList[i], CL_FALSE, 0,
                                 this->sizeofType * this->elements, hostData, 0, nullptr, nullptr);
    } else {
      err = clEnqueueWriteBuffer(this->manager->command_queue[i], this->bufferList[i], CL_FALSE,
                                 this->sizeofType * offsets[i], this->sizeofType * this->elements,
                                 hostData, 0, nullptr, nullptr);
    }

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to enqueue write "
                     "buffer command! Error code: " << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  clWaitForEvents((cl_uint) this->manager->num_devices, actionDone);

  for (size_t i = 0; i < this->manager->num_devices; i++) {
    clReleaseEvent(actionDone[i]);
  }
}

void OCLClonedBuffer::readFromBuffer(void* hostData, size_t* offsets, size_t* ranges) {
  cl_int err;
  cl_event* actionDone = new cl_event[this->manager->num_devices];

  // read data back
  for (size_t i = 0; i < this->manager->num_devices; i++) {
    if (offsets == nullptr) {
      err = clEnqueueReadBuffer(this->manager->command_queue[i], this->bufferList[i], CL_FALSE, 0,
                                this->sizeofType * this->elements, hostData, 0, nullptr,
                                &(actionDone[i]));
    } else {
      err = clEnqueueReadBuffer(this->manager->command_queue[i], this->bufferList[i], CL_FALSE,
                                this->sizeofType * offsets[i], this->sizeofType * ranges[i],
                                static_cast<char*>(hostData) + (this->sizeofType * offsets[i]), 0,
                                nullptr, &(actionDone[i]));
    }

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to enqueue read buffer command! "
                     "Error code: " << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  clWaitForEvents((cl_uint) this->manager->num_devices, actionDone);

  for (size_t i = 0; i < this->manager->num_devices; i++) {
    clReleaseEvent(actionDone[i]);
  }
}

// read/write-flags are missing
void OCLClonedBuffer::initializeBuffer(void* initialValues, size_t sizeofType, size_t elements) {
  cl_int err;
  cl_mem* bufferList = new cl_mem[manager->num_devices];

  for (size_t i = 0; i < manager->num_devices; i++) {
    if (initialValues != nullptr) {
      bufferList[i] = clCreateBuffer(manager->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                     sizeofType * elements, initialValues, &err);
    } else {
      bufferList[i] =
          clCreateBuffer(manager->context, CL_MEM_READ_ONLY, sizeofType * elements, nullptr, &err);
    }

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Could not allocate buffer! "
                     "Error code: " << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  this->bufferList = bufferList;
  this->sizeofType = sizeofType;
  this->elements = elements;
  this->initialized = true;
}

void OCLClonedBuffer::freeBuffer() {
  if (!this->initialized) {
    return;
  }

  if (this->bufferList == nullptr) {
    std::stringstream errorString;
    errorString << "OCL Error: OCLClonedBuffer in partially initialized state: "
                   "buffer list is null" << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  for (size_t i = 0; i < this->manager->num_devices; i++) {
    if (this->bufferList[i] != nullptr) {
      clReleaseMemObject(this->bufferList[i]);
      this->bufferList[i] = nullptr;
    } else {
      std::stringstream errorString;
      errorString << "OCL Error: OCLClonedBuffer in partially initialized state: "
                     "device buffer is null" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  delete[] this->bufferList;
  this->bufferList = nullptr;
  this->initialized = false;
}

}  // namespace base
}  // namespace sgpp

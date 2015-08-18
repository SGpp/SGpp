/*
 * OCLMemory.cpp
 *
 *  Created on: Mar 30, 2015
 *      Author: pfandedd
 */

#include "../../../../../datadriven/src/sgpp/datadriven/opencl/OCLClonedBuffer.hpp"

#include <sstream>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/exception/operation_exception.hpp>

namespace SGPP {
  namespace base {

    OCLClonedBuffer::OCLClonedBuffer(std::shared_ptr<OCLManager> manager) :
      manager(manager) {
      initialized = false;
      bufferList = nullptr;
      sizeofType = 0;
      elements = 0;
    }

    bool OCLClonedBuffer::isInitialized() {
      return this->initialized;
    }

    cl_mem* OCLClonedBuffer::getBuffer(size_t deviceNumber) {
      return &(this->bufferList[deviceNumber]);
    }

    //TODO: current multidevice strategy: allocate everything everywere, use only range specified for device

    void OCLClonedBuffer::writeToBuffer(void* hostData, size_t* offsets) {
      cl_int err;
      cl_event* actionDone = new cl_event[this->manager->num_devices];

      //TODO: fix bug, not waiting for event
      for (size_t i = 0; i < this->manager->num_devices; i++) {
        if (offsets == nullptr) {
          err = clEnqueueWriteBuffer(this->manager->command_queue[i], this->bufferList[i],
                                     CL_FALSE, 0, this->sizeofType * this->elements, hostData, 0, nullptr, nullptr);
        } else {
          err = clEnqueueWriteBuffer(this->manager->command_queue[i], this->bufferList[i],
                                     CL_FALSE, this->sizeofType * offsets[i], this->sizeofType * this->elements, hostData, 0, nullptr, nullptr);
        }

        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to enqueue write buffer command! Error code: " << err << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
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
          err = clEnqueueReadBuffer(this->manager->command_queue[i], this->bufferList[i],
                                    CL_FALSE, 0, this->sizeofType * this->elements, hostData, 0, nullptr, &(actionDone[i]));
        } else {
          err = clEnqueueReadBuffer(this->manager->command_queue[i], this->bufferList[i],
                                    CL_FALSE, this->sizeofType * offsets[i], this->sizeofType * ranges[i],
                                    static_cast<char*>(hostData) + (this->sizeofType * offsets[i]), 0, nullptr, &(actionDone[i]));
        }

        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to enqueue read buffer command! Error code: " << err << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }
      }

      clWaitForEvents((cl_uint) this->manager->num_devices, actionDone);

      for (size_t i = 0; i < this->manager->num_devices; i++) {
        clReleaseEvent(actionDone[i]);
      }
    }

    //read/write-flags are missing
    void OCLClonedBuffer::initializeBuffer(void* initialValues, size_t sizeofType, size_t elements) {
      cl_int err;
      cl_mem* bufferList = new cl_mem[manager->num_devices];

      for (size_t i = 0; i < manager->num_devices; i++) {
        if (initialValues != nullptr) {
          bufferList[i] = clCreateBuffer(manager->context,
                                         CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeofType * elements, initialValues, &err);
        } else {
          bufferList[i] = clCreateBuffer(manager->context,
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
      this->initialized = true;
    }

    void OCLClonedBuffer::freeBuffer() {
      if (!this->initialized) {
        return;
      }

      for (size_t i = 0; i < this->manager->num_devices; i++) {
        if (this->bufferList[i] != nullptr) {
          clReleaseMemObject(this->bufferList[i]);
          this->bufferList[i] = nullptr;
        }
      }

      delete[] this->bufferList;
      this->initialized = false;
    }

  }
}


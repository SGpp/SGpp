/*
 * OCLMemory.hpp
 *
 *  Created on: Mar 27, 2015
 *      Author: pfandedd
 */

#pragma once

#include <CL/cl.h>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>

namespace SGPP {
namespace base {

template<typename T>
class OCLClonedBufferSD {
private:
    std::shared_ptr<OCLDevice> device;
    bool initialized;
    cl_mem buffer;
//    size_t sizeofType;
    size_t elements;
    std::vector<T> hostData;

public:

    OCLClonedBufferSD(std::shared_ptr<base::OCLDevice> device) :
            device(device), initialized(false), buffer(nullptr), elements(0) {
    }

    ~OCLClonedBufferSD() {
        this->freeBuffer();
    }

    bool isInitialized() {
        return this->initialized;
    }

    cl_mem *getBuffer() {
        if (!this->initialized) {
            std::stringstream errorString;
            errorString << "OCL Error: Buffer not initialized: " << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }
        return &this->buffer;
    }

    size_t size() {
        if (!this->initialized) {
            std::stringstream errorString;
            errorString << "OCL Error: Buffer not initialized: " << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }
        return this->elements;
    }

    std::vector<T> &getHostPointer() {
        if (!this->initialized) {
            std::stringstream errorString;
            errorString << "OCL Error: Buffer not initialized: " << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }
        return hostData;
    }

    void writeToBuffer() {
        if (!this->initialized) {
            std::stringstream errorString;
            errorString << "OCL Error: Buffer not initialized: " << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }

        cl_int err;

        err = clEnqueueWriteBuffer(device->commandQueue, this->buffer,
        CL_TRUE, 0, sizeof(T) * this->elements, hostData.data(), 0, nullptr, nullptr);

        if (err != CL_SUCCESS) {
            std::stringstream errorString;
            errorString << "OCL Error: Failed to enqueue write buffer command! Error code: " << err << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }
    }

    void readFromBuffer() {
        if (!this->initialized) {
            std::stringstream errorString;
            errorString << "OCL Error: Buffer not initialized: " << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }

        cl_int err;

        std::cout << "cc in buffer: " << device->commandQueue << std::endl;
        err = clEnqueueReadBuffer(device->commandQueue, this->buffer,
        CL_TRUE, 0, sizeof(T) * this->elements,
                static_cast<void *>(hostData.data()), 0, nullptr, nullptr);

        if (err != CL_SUCCESS) {
            std::stringstream errorString;
            errorString << "OCL Error: Failed to enqueue read buffer command! Error code: " << err << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }
    }

    //TODO: might need to set the correct flags
    void initializeBuffer(size_t elements, std::vector<T> initialValues = std::vector<T>()) {
        cl_int err;

        hostData.resize(elements);

        if (initialValues.size() != 0) {
            if (initialValues.size() != elements) {
                std::stringstream errorString;
                errorString << "OCL Error: Size of initial values vector does not match specified size! Error code: "
                        << err << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
            }
            this->buffer = clCreateBuffer(device->context,
            CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(T) * elements, initialValues.data(), &err);
        } else {
            this->buffer = clCreateBuffer(device->context,
            CL_MEM_READ_WRITE, sizeof(T) * elements, nullptr, &err);
        }

        if (err != CL_SUCCESS) {
            std::stringstream errorString;
            errorString << "OCL Error: Could not allocate buffer! Error code: " << err << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }

        this->elements = elements;
        this->initialized = true;
    }

    void freeBuffer() {
        if (!this->initialized) {
            //permitted to allow for easy kernel resets
            return;
        }

        if (buffer != nullptr) {
            clReleaseMemObject(buffer);
            buffer = nullptr;
        } else {
            std::stringstream errorString;
            errorString << "OCL Error: could not free OCLClonedBufferSD" << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }
        this->initialized = false;
        this->elements = 0;
    }

};

}
}


/*
 * OCLReadOnlyBuffer.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: leiterrl
 */

#include "../../../../../datadriven/src/sgpp/datadriven/opencl/OCLReadOnlyBuffer.hpp"

#include <sstream>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/exception/operation_exception.hpp>
#include <cstring>

namespace SGPP
{
    namespace base
    {

        OCLReadOnlyBuffer::OCLReadOnlyBuffer(std::shared_ptr<OCLManager> manager) : m_manager(manager)
        {

            m_initialized = false;
            m_bufferList = nullptr;
            m_mappedHostBuffer = nullptr;
            m_sizeofType = 0;
            m_elements = 0;

        }

        OCLReadOnlyBuffer::~OCLReadOnlyBuffer()
        {
            this->freeBuffer();
        }

        bool OCLReadOnlyBuffer::isInitialized()
        {
            return this->m_initialized;
        }

        cl_mem* OCLReadOnlyBuffer::getBuffer(size_t deviceNumber)
        {
            return &(m_hostBuffer);
        }

        //TODO: current multidevice strategy: allocate everything everywere, use only range specified for device

        void OCLReadOnlyBuffer::writeToBuffer(void* hostData)
        {
            cl_int err;

            m_mappedHostBuffer = clEnqueueMapBuffer(m_manager->command_queue[0], m_hostBuffer, CL_TRUE,
                    CL_MAP_WRITE, 0, m_sizeofType * m_elements, 0, nullptr, nullptr, &err);

            if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to enqueue map buffer command when trying to write! Error code: " << err << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
            }

            memcpy(m_mappedHostBuffer, hostData, m_sizeofType*m_elements);

            err = clEnqueueUnmapMemObject(m_manager->command_queue[0], m_hostBuffer, m_mappedHostBuffer, 0, nullptr, nullptr);

            if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to enqueue unmap memobject command when trying to write! Error code: " << err << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
            }

        }

        void OCLReadOnlyBuffer::readFromBuffer(void* hostData)
        {
            cl_int err;

            m_mappedHostBuffer = clEnqueueMapBuffer(m_manager->command_queue[0], m_hostBuffer, CL_TRUE,
                    CL_MAP_READ, 0, m_sizeofType * m_elements, 0, nullptr, nullptr, &err);

            if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to enqueue map buffer command when trying to read! Error code: " << err << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
            }

            memcpy(hostData, m_mappedHostBuffer, m_sizeofType*m_elements);

            err = clEnqueueUnmapMemObject(m_manager->command_queue[0], m_hostBuffer, m_mappedHostBuffer, 0, nullptr, nullptr);

            if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to enqueue unmap memobject command when trying to read! Error code: " << err << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
            }
        }

        //read/write-flags are missing
        void OCLReadOnlyBuffer::initializeBuffer(void* initialValues, size_t sizeofType, size_t elements)
        {
            cl_int err;
            cl_mem* bufferList = new cl_mem[m_manager->num_devices];

            m_hostBuffer = clCreateBuffer(m_manager->context,
                CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeofType * elements, nullptr, &err);

            if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Could not allocate buffer! Error code: " << err << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
            }

            this->m_bufferList = bufferList;
            this->m_sizeofType = sizeofType;
            this->m_elements = elements;
            this->m_initialized = true;

            if (initialValues != nullptr)
            {
                writeToBuffer(initialValues);
            }
        }

        void OCLReadOnlyBuffer::freeBuffer()
        {
            if (!this->m_initialized) {
                return;
            }

            if (this->m_bufferList == nullptr) {
                std::stringstream errorString;
                errorString << "OCL Error: OCLClonedBuffer in partially m_initialized state: buffer list is null" << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
            }

            for (size_t i = 0; i < this->m_manager->num_devices; i++) {
                if (this->m_bufferList[i] != nullptr) {
                    clReleaseMemObject(this->m_bufferList[i]);
                    this->m_bufferList[i] = nullptr;
                } else {
                    std::stringstream errorString;
                    errorString << "OCL Error: OCLClonedBuffer in partially m_initialized state: device buffer is null"
                            << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }
            }

            delete[] this->m_bufferList;
            this->m_bufferList = nullptr;
            this->m_initialized = false;
        }

    }
}

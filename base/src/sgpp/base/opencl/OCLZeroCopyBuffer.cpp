// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sstream>
#include <cstring>

#include <sgpp/base/opencl/OCLZeroCopyBuffer.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

namespace sgpp {
namespace base {

OCLZeroCopyBuffer::OCLZeroCopyBuffer(std::shared_ptr<OCLManager> manager)
    : m_manager(manager), m_hostBuffer(nullptr) {
  m_initialized = false;
  m_mappedHostBuffer = nullptr;
  m_sizeofType = 0;
  m_elements = 0;
  m_readOnly = false;
}

OCLZeroCopyBuffer::~OCLZeroCopyBuffer() { this->freeBuffer(); }

bool OCLZeroCopyBuffer::isInitialized() { return this->m_initialized; }

cl_mem* OCLZeroCopyBuffer::getBuffer(size_t deviceNumber) { return &(m_hostBuffer); }

void OCLZeroCopyBuffer::writeToBuffer(void* hostData) {
  cl_int err;

  m_mappedHostBuffer =
      clEnqueueMapBuffer(m_manager->command_queue[0], m_hostBuffer, CL_TRUE, CL_MAP_WRITE, 0,
                         m_sizeofType * m_elements, 0, nullptr, nullptr, &err);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Failed to enqueue map buffer command when "
                   "trying to write! Error code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  if (m_mappedHostBuffer == nullptr)
    throw std::runtime_error("OCLZeroCopyBuffer::writeToBuffer mappedHostBuffer == NULL");

  /*for ( int i = 0; i < m_elements; i++)
  {
      printf("hostData[%i]: %i \n", i, *((uint32_t*)hostData + 0x4*i));
  }*/

  // uint32_t* testArray = new uint32_t[m_elements];
  // memcpy(m_mappedHostBuffer, testArray, m_sizeofType*m_elements);

  memcpy(m_mappedHostBuffer, hostData, m_sizeofType * m_elements);

  err = clEnqueueUnmapMemObject(m_manager->command_queue[0], m_hostBuffer, m_mappedHostBuffer, 0,
                                nullptr, nullptr);
  m_mappedHostBuffer = nullptr;

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Failed to enqueue unmap memobject command "
                   "when trying to write! Error code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }
}

void OCLZeroCopyBuffer::readFromBuffer(void* hostData) {
  cl_int err;

  m_mappedHostBuffer =
      clEnqueueMapBuffer(m_manager->command_queue[0], m_hostBuffer, CL_TRUE, CL_MAP_READ, 0,
                         m_sizeofType * m_elements, 0, nullptr, nullptr, &err);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Failed to enqueue map buffer command "
                   "when trying to read! Error code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  if (m_mappedHostBuffer == nullptr)
    throw std::runtime_error(
        "OCLZeroCopyBuffer::writeToBuffer "
        "mappedHostBuffer == NULL");

  memcpy(hostData, m_mappedHostBuffer, m_sizeofType * m_elements);

  err = clEnqueueUnmapMemObject(m_manager->command_queue[0], m_hostBuffer, m_mappedHostBuffer, 0,
                                nullptr, nullptr);
  m_mappedHostBuffer = nullptr;

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Failed to enqueue unmap memobject command "
                   "when trying to read! Error code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }
}

// read/write-flags are missing
void OCLZeroCopyBuffer::initializeBuffer(void* initialValues, size_t sizeofType, size_t elements,
                                         bool readOnly) {
  cl_int err;

  cl_mem_flags flags = CL_MEM_ALLOC_HOST_PTR;

  if (readOnly) {
    flags |= CL_MEM_READ_ONLY;
  } else {
    flags |= CL_MEM_READ_WRITE;
  }

  m_hostBuffer = clCreateBuffer(m_manager->context, flags, sizeofType * elements, nullptr, &err);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Could not allocate buffer! "
                   "Error code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  this->m_sizeofType = sizeofType;
  this->m_elements = elements;
  this->m_initialized = true;
  this->m_readOnly = readOnly;

  if (initialValues != nullptr) {
    writeToBuffer(initialValues);
  }
}

void OCLZeroCopyBuffer::freeBuffer() {
  if (!this->m_initialized) {
    return;
  }

  if (m_hostBuffer != nullptr) {
    clReleaseMemObject(m_hostBuffer);
    m_hostBuffer = nullptr;
  }

  this->m_initialized = false;
}

}  // namespace base
}  // namespace sgpp

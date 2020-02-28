// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <CL/cl.h>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>

namespace sgpp {
namespace base {

class OCLZeroCopyBuffer {
 private:
  std::shared_ptr<OCLManager> m_manager;
  bool m_initialized;
  size_t m_sizeofType;
  size_t m_elements;
  cl_mem m_hostBuffer;
  void* m_mappedHostBuffer;
  bool m_readOnly;

 public:
  explicit OCLZeroCopyBuffer(std::shared_ptr<OCLManager> manager);
  ~OCLZeroCopyBuffer();

  bool isInitialized();

  cl_mem* getBuffer(size_t deviceNumber);

  void* getMappedBuffer();

  void writeToBuffer(void* hostData);

  void readFromBuffer(void* hostData);

  void initializeBuffer(void* initialValues, size_t sizeofType, size_t elements,
                        bool readOnly);

  void freeBuffer();
};

}  // namespace base
}  // namespace sgpp

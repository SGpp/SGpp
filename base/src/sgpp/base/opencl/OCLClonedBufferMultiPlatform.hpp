/*
 * OCLMemory.hpp
 *
 *  Created on: Mar 27, 2015
 *      Author: pfandedd
 */

#pragma once

#include <CL/cl.h>

#include <sgpp/globaldef.hpp>

#include "OCLManagerMultiPlatform.hpp"

namespace SGPP {
  namespace base {

    class OCLClonedBufferMultiPlatform {
      public:
        OCLManagerMultiPlatform& manager;
        bool initialized;
        std::map<cl_platform_id, cl_mem *> platformBufferList;
        size_t sizeofType;
        size_t elements;

      public:

        OCLClonedBufferMultiPlatform(OCLManagerMultiPlatform& manager);

        bool isInitialized();

        cl_mem* getBuffer(cl_platform_id platformId, size_t deviceIndex);

        void writeToBuffer(void* hostData, size_t* offsets = nullptr);

        void readFromBuffer(void* hostData, size_t* offsets = nullptr, size_t* ranges = nullptr);

        void initializeBuffer(void* initialValues, size_t sizeofType, size_t elements);

        void freeBuffer();

    };

  }
}



/*
 * OCLReadOnlyBuffer.h
 *
 *  Created on: Oct 9, 2015
 *      Author: leiterrl
 */

#pragma once

#include <CL/cl.h>

#include <sgpp/globaldef.hpp>

#include "../../../../../datadriven/src/sgpp/datadriven/opencl/OCLManager.hpp"

namespace SGPP
{
    namespace base
    {

        class OCLReadOnlyBuffer
        {
        private:
            std::shared_ptr<OCLManager> m_manager;
            bool m_initialized;
            cl_mem* m_bufferList;
            size_t m_sizeofType;
            size_t m_elements;
            cl_mem m_hostBuffer;
            void* m_mappedHostBuffer;

        public:
            OCLReadOnlyBuffer(std::shared_ptr<OCLManager> manager);
            ~OCLReadOnlyBuffer();

            bool isInitialized();

            cl_mem* getBuffer(size_t deviceNumber);

            void* getMappedBuffer();

            void writeToBuffer(void* hostData);

            void readFromBuffer(void* hostData);

            void initializeBuffer(void* initialValues, size_t sizeofType, size_t elements);

            void freeBuffer();
        };

    } /* namespace base */
} /* namespace SGPP */

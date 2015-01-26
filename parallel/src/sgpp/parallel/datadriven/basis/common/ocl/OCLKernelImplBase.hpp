/* ****************************************************************************
* Copyright (C) 2010-2013 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef OCLKERNELIMPLBASE_HPP
#define OCLKERNELIMPLBASE_HPP

#ifdef USEOCL
// include OpenCL
#include <CL/cl.h>

#include <string.h>

#include <sgpp/parallel/datadriven/basis/common/KernelBase.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>
#include <limits>

#ifdef USEOCL_NVIDIA
#define USEOCL_LOCAL_MEMORY
#endif

#ifdef USEOCL_AMD
#ifndef USEOCL_CPU
#define USEOCL_LOCAL_MEMORY
#endif
#endif

namespace sg {
  namespace parallel {
    class OCLKernelImplBase {
      public:
        void resetKernel();
        static inline size_t getChunkGridPoints() {
          return ocl_local_size;
        }

        static inline size_t getChunkDataPoints() {
          return ocl_local_size;
        }
      protected:
        OCLKernelImplBase();

        ~OCLKernelImplBase();
        void releaseGridBuffers();

        void releaseDataBuffers();

        void releaseKernelsAndPrograms();


        cl_int err;
        cl_platform_id platform_id;
        cl_platform_id* platform_ids;
        cl_device_id* device_ids;
        cl_uint num_platforms;
        cl_uint num_devices;
        cl_context context;

        cl_command_queue* command_queue;

        cl_mem* clData;
        cl_mem* clLevel;
        cl_mem* clIndex;
        cl_mem* clMask;
        cl_mem* clOffset;

        // use pinned memory (on host and device) to speed up data transfers from/to GPU
        cl_mem* clDevGrid;
        cl_mem* clDevTmp;
        cl_mem clPinnedGrid;
        cl_mem clPinnedTmp;

        cl_kernel* kernel_multTrans;
        cl_kernel* kernel_mult;

        static unsigned int ocl_local_size;
    };
  }
}

#endif // USEOCL

#endif // OCLKERNELIMPLBASE_HPP

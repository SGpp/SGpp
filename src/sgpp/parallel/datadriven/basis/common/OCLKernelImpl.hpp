/* ****************************************************************************
* Copyright (C) 2010-2013 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef OCLKERNELIMPL_HPP
#define OCLKERNELIMPL_HPP

// include OpenCL
#include <CL/cl.h>

#include <string.h>

#include "parallel/datadriven/basis/common/KernelBase.hpp"
#include "parallel/tools/PartitioningTool.hpp"

#define MAX_OCL_DEVICE_COUNT 8
#ifndef USEOCL_LOCAL_WORKGROUP_SIZE
#define OCL_SGPP_LOCAL_WORKGROUP_SIZE 64
#else
#define OCL_SGPP_LOCAL_WORKGROUP_SIZE USEOCL_LOCAL_WORKGROUP_SIZE
#endif

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

    template<typename OCLBasisType>
    class OCLKernelImpl {
      protected:
        cl_int err;
        cl_platform_id platform_id;
        cl_platform_id* platform_ids;
        cl_device_id* device_ids;
        cl_uint num_platforms;
        cl_uint num_devices;
        cl_context context;

        cl_command_queue command_queue[MAX_OCL_DEVICE_COUNT];

        cl_mem clData[MAX_OCL_DEVICE_COUNT];
        cl_mem clLevel[MAX_OCL_DEVICE_COUNT];
        cl_mem clIndex[MAX_OCL_DEVICE_COUNT];
        cl_mem clMask[MAX_OCL_DEVICE_COUNT];
        cl_mem clOffset[MAX_OCL_DEVICE_COUNT];

        cl_kernel kernel_multTrans[MAX_OCL_DEVICE_COUNT];
        cl_kernel kernel_mult[MAX_OCL_DEVICE_COUNT];

      public:
        OCLKernelImpl() {
          // initialize arrays
          for (int i = 0; i < MAX_OCL_DEVICE_COUNT; i++) {
            command_queue[i] = NULL;

            clData[i] = NULL;
            clLevel[i] = NULL;
            clIndex[i] = NULL;
            clMask[i] = NULL;
            clOffset[i] = NULL;

            kernel_mult[i] = NULL;
            kernel_multTrans[i] = NULL;
          }

          // determine number of available OpenCL platforms
          err = clGetPlatformIDs(0, NULL, &num_platforms);

          if (err != CL_SUCCESS) {
            std::cout << "OCL Error: Unable to get number of OpenCL platforms. Error Code: " << err << std::endl;
          }

          std::cout << "OCL Info: " << num_platforms << " OpenCL Platforms have been found" << std::endl;

          // get available platforms
          platform_ids = new cl_platform_id[num_platforms];
          err = clGetPlatformIDs(num_platforms, platform_ids, NULL);

          if (err != CL_SUCCESS) {
            std::cout << "OCL Error: Unable to get Platform ID. Error Code: " << err << std::endl;
          }

          for (cl_uint ui = 0; ui < num_platforms; ui++) {
            char vendor_name[128] = {0};
            err = clGetPlatformInfo(platform_ids[ui], CL_PLATFORM_VENDOR, 128 * sizeof(char), vendor_name, NULL);

            if (CL_SUCCESS != err) {
              std::cout << "OCL Error: Can't get platform vendor!" << std::endl;
            } else {
              if (vendor_name != NULL) {
                std::cout << "OCL Info: Platform " << ui << " vendor name: " << vendor_name << std::endl;
              }

#ifdef USEOCL_INTEL

              if (strcmp(vendor_name, "Intel(R) Corporation") == 0) {
#ifdef USEOCL_CPU
                std::cout << "OCL Info: Using CPU Platform: " << vendor_name << std::endl;
#else
                std::cout << "OCL Info: Using GPU Platform: " << vendor_name << std::endl;
#endif
                platform_id = platform_ids[ui];
              }

#endif
#ifdef USEOCL_AMD

              if (strcmp(vendor_name, "Advanced Micro Devices, Inc.") == 0) {
#ifdef USEOCL_CPU
                std::cout << "OCL Info: Using CPU Platform: " << vendor_name << std::endl;
#else
                std::cout << "OCL Info: Using GPU Platform: " << vendor_name << std::endl;
#endif
                platform_id = platform_ids[ui];
              }

#endif
#ifdef USEOCL_NVIDIA

              if (strcmp(vendor_name, "NVIDIA Corporation") == 0) {
                std::cout << "OCL Info: Using GPU Platform: " << vendor_name << std::endl;
                platform_id = platform_ids[ui];
              }

#endif
            }
          }

          std::cout << std::endl;

          // Find out how many devices there are
#ifdef USEOCL_INTEL
          device_ids = new cl_device_id[1];
#ifdef USEOCL_CPU
          err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, 1, device_ids, NULL);
#else
          err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, device_ids, NULL);
#endif

          if (err != CL_SUCCESS) {
            std::cout << "OCL Error: Unable to get Device ID. Error Code: " << err << std::endl;
          }

          num_devices = 1;
#endif
#ifdef USEOCL_AMD
          device_ids = new cl_device_id[MAX_OCL_DEVICE_COUNT];
          err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, MAX_OCL_DEVICE_COUNT, device_ids, &num_devices);

          if (num_devices == 0) {
            std::cout << "OCL Error: NO GPU OpenCL devices have been found!" << std::endl;
          }

          // set max number of devices
          if (num_devices > MAX_OCL_DEVICE_COUNT) {
            num_devices = MAX_OCL_DEVICE_COUNT;
          }

#ifdef USEOCL_CPU
          err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, 2, device_ids, NULL);
#else
          err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 2, device_ids, NULL);
#endif

          if (err != CL_SUCCESS) {
            std::cout << "OCL Error: Unable to get Device ID. Error Code: " << err << std::endl;
          }

#endif
#ifdef USEOCL_NVIDIA
          err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, MAX_OCL_DEVICE_COUNT, NULL, &num_devices);

          if (num_devices == 0) {
            std::cout << "OCL Error: NO GPU OpenCL devices have been found!" << std::endl;
          }

          // set max number of devices
          if (num_devices > MAX_OCL_DEVICE_COUNT) {
            num_devices = MAX_OCL_DEVICE_COUNT;
          }

          device_ids = new cl_device_id[num_devices];
          err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, num_devices, device_ids, NULL);

          if (err != CL_SUCCESS) {
            std::cout << "OCL Error: Unable to get Device ID. Error Code: " << err << std::endl;
          }

#endif
          std::cout << "OCL Info: " << num_devices << " OpenCL devices have been found!" << std::endl;

          // Create OpenCL context
          context = clCreateContext(0, num_devices, device_ids, NULL, NULL, &err);

          if (err != CL_SUCCESS) {
            std::cout << "OCL Error: Failed to create OpenCL context! Error Code: " << err << std::endl;
          }

          // Creating the command queues
          for (size_t i = 0; i < num_devices; i++) {
            // TODO FIXME whats the difference here?
#ifdef USEOCL_CPU
            command_queue[i] = clCreateCommandQueue(context, device_ids[i], CL_QUEUE_PROFILING_ENABLE, &err);
#else
            command_queue[i] = clCreateCommandQueue(context, device_ids[i], CL_QUEUE_PROFILING_ENABLE, &err);
#endif

            if (err != CL_SUCCESS) {
              std::cout << "OCL Error: Failed to create command queue! Error Code: " << err << std::endl;
            }
          }

          std::cout << "OCL Info: Successfully initialized OpenCL (local workgroup size: " << OCL_SGPP_LOCAL_WORKGROUP_SIZE << ")" << std::endl << std::endl;
        }

        ~OCLKernelImpl() {
          releaseDataBuffers();
          releaseGridBuffers();
          releaseKernelsAndPrograms();

          // release command queue
          for (size_t i = 0; i < num_devices; i++) {
            if (command_queue[i]) {
              clReleaseCommandQueue(command_queue[i]);
            }
          }

          // release context
          clReleaseContext(context);

          delete[] device_ids;
        }

        inline void initOCLBuffers(
          sg::base::DataMatrix* level,
          sg::base::DataMatrix* index,
          sg::base::DataMatrix* mask,
          sg::base::DataMatrix* offset,
          sg::base::DataMatrix* dataset) {
          size_t storageSize = level->getSize();

          if (level != NULL && clLevel[0] == NULL) {
            for (size_t i = 0; i < num_devices; i++) {
              clLevel[i] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * storageSize, level->getPointer(), NULL);
            }
          }

          if (index != NULL && clIndex[0] == NULL) {
            for (size_t i = 0; i < num_devices; i++) {
              clIndex[i] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * storageSize, index->getPointer(), NULL);
            }
          }

          if (mask != NULL && clMask[0] == NULL) {
            for (size_t i = 0; i < num_devices; i++) {
              clMask[i] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * storageSize, mask->getPointer(), NULL);
            }
          }

          if (offset != NULL && clOffset[0] == NULL) {
            for (size_t i = 0; i < num_devices; i++) {
              clOffset[i] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * storageSize, offset->getPointer(), NULL);
            }
          }

          if (clData[0] == NULL) { // use first element as indicator if data has been already copied to device
            for (size_t i = 0; i < num_devices; i++) {
              clData[i] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * dataset->getSize(), dataset->getPointer(), NULL);
            }
          }
        }

        double multImpl(
          sg::base::DataMatrix* level,
          sg::base::DataMatrix* index,
          sg::base::DataMatrix* mask,
          sg::base::DataMatrix* offset,
          sg::base::DataMatrix* dataset,
          sg::base::DataVector& alpha,
          sg::base::DataVector& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {
          // check if there is something to do at all
          if (!(end_index_grid > start_index_grid && end_index_data > start_index_data)) {
            return 0.0;
          }

          double* ptrResult = result.getPointer();
          double time = 0.0;

          if (kernel_mult[0] == NULL) {
            OCLBasisType basis;
            size_t dims = dataset->getNrows();
            basis.createMult(dims, context, num_devices, device_ids, kernel_mult);
          }

          initOCLBuffers(level, index, mask, offset, dataset);

          cl_mem clAlpha[MAX_OCL_DEVICE_COUNT], clResult[MAX_OCL_DEVICE_COUNT];

          for (size_t i = 0; i < num_devices; i++) {
            clAlpha[i] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * alpha.getSize(), alpha.getPointer(), NULL);
            clResult[i] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(double) * result.getSize(), result.getPointer(), NULL);
          }

          // determine best fit
          size_t gpu_start_index_data[MAX_OCL_DEVICE_COUNT];
          size_t gpu_end_index_data[MAX_OCL_DEVICE_COUNT];

          for (size_t gpu_num = 0; gpu_num < num_devices; gpu_num++) {
            sg::parallel::PartitioningTool::getPartitionSegment(start_index_data, end_index_data, num_devices, gpu_num, &gpu_start_index_data[gpu_num], &gpu_end_index_data[gpu_num], OCL_SGPP_LOCAL_WORKGROUP_SIZE);
          }

          // set kernel arguments
          cl_uint clResultSize = (cl_uint)(result.getSize());
          cl_uint gpu_start_grid = start_index_grid;
          cl_uint gpu_end_grid = end_index_grid;

          for (size_t i = 0; i < num_devices; i++) {
            cl_uint gpu_start_data = (cl_uint)gpu_start_index_data[i];
            cl_uint gpu_end_data = (cl_uint)gpu_end_index_data[i];

            if (gpu_end_data > gpu_start_data) {
              if (clSetKernelArg(kernel_mult[i], 0, sizeof(cl_mem), &clLevel[i]) ||
                  clSetKernelArg(kernel_mult[i], 1, sizeof(cl_mem), &clIndex[i]) ||
                  clSetKernelArg(kernel_mult[i], 2, sizeof(cl_mem), &clMask[i]) || // only needed for masked version of modlinear
                  clSetKernelArg(kernel_mult[i], 3, sizeof(cl_mem), &clOffset[i]) || //only needed for masked version of modlinear
                  clSetKernelArg(kernel_mult[i], 4, sizeof(cl_mem), &clData[i]) ||
                  clSetKernelArg(kernel_mult[i], 5, sizeof(cl_mem), &clAlpha[i]) ||
                  clSetKernelArg(kernel_mult[i], 6, sizeof(cl_mem), &clResult[i]) ||
                  clSetKernelArg(kernel_mult[i], 7, sizeof(cl_uint), &clResultSize) || // resultsize == number of entries in dataset
                  clSetKernelArg(kernel_mult[i], 8, sizeof(cl_uint), &gpu_start_grid) ||
                  clSetKernelArg(kernel_mult[i], 9, sizeof(cl_uint), &gpu_end_grid) != CL_SUCCESS) {
                std::cout << "OCL Error: Failed to create kernel Args for mult!" << std::endl;
                return 0.0;
              }
            }
          }

          cl_event clTimings[MAX_OCL_DEVICE_COUNT];
          cl_event GPUDone[MAX_OCL_DEVICE_COUNT];

          // enqueue kernel
          size_t local = OCL_SGPP_LOCAL_WORKGROUP_SIZE;
          size_t active_devices = 0;

          for (size_t i = 0; i < num_devices; i++) {
            size_t rangeSize = gpu_end_index_data[i] - gpu_start_index_data[i];

            if (rangeSize > 0) {
              err = clEnqueueNDRangeKernel(command_queue[i], kernel_mult[i], 1, &gpu_start_index_data[i], &rangeSize, &local, 0, NULL, &(clTimings[i]));

              if (active_devices != i) {
                std::cout << "OCL Error: Splitting up calculations for multiple GPUs is erroneous, only the last chunks may handle 0 entries. syncing will be not correct.";
              }

              active_devices++;

              if (err != CL_SUCCESS) {
                std::cout << "OCL Error: Failed to enqueue kernel command! Error Code: " << err << std::endl;
                return 0.0;
              }
            }
          }

          // read data back
          for (size_t i = 0; i < num_devices; i++) {
            size_t rangeSize = gpu_end_index_data[i] - gpu_start_index_data[i];

            if (rangeSize > 0) {
              size_t offset = gpu_start_index_data[i];
              err = clEnqueueReadBuffer(command_queue[i], clResult[i], CL_FALSE, sizeof(double) * offset, sizeof(double) * rangeSize, &(ptrResult[offset]), 0, NULL, &(GPUDone[i]));

              if (err != CL_SUCCESS) {
                std::cout << "OCL Error: Failed to enqueue read buffer command (mult)! Error Code: " << err << std::endl;
                return 0.0;
              }
            }
          }

          // sync GPUs
          clWaitForEvents((cl_uint)active_devices, GPUDone);

          // determine kernel execution time
          for (size_t i = 0; i < num_devices; i++) {
            double tmpTime;
            cl_ulong startTime, endTime;
            startTime = endTime = 0;

            if (gpu_end_index_data[i] > gpu_start_index_data[i]) {
              err = clGetEventProfilingInfo(clTimings[i], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, NULL);

              if (err != CL_SUCCESS) {
                std::cout << "OCL Error: Failed to read start-time from command queue! Error Code: " << err << std::endl;
              }

              err = clGetEventProfilingInfo(clTimings[i], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, NULL);

              if (err != CL_SUCCESS) {
                std::cout << "OCL Error: Failed to read end-time from command queue! Error Code: " << err << std::endl;
              }
            }

            tmpTime = (double)(endTime - startTime);
            tmpTime *= 1e-9;

            if (tmpTime > time) {
              time = tmpTime;
            }
          }

          // clean up
          for (size_t i = 0; i < num_devices; i++) {
            clReleaseMemObject(clAlpha[i]);
            clReleaseMemObject(clResult[i]);

            if (gpu_end_index_data[i] > gpu_start_index_data[i]) {
              clReleaseEvent(clTimings[i]);
              clReleaseEvent(GPUDone[i]);
            }
          }

          return time;
        }

        double multTransposeImpl(
          sg::base::DataMatrix* level,
          sg::base::DataMatrix* index,
          sg::base::DataMatrix* mask,
          sg::base::DataMatrix* offset,
          sg::base::DataMatrix* dataset,
          sg::base::DataVector& source,
          sg::base::DataVector& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {
          // check if there is something to do at all
          if ( !(end_index_grid > start_index_grid && end_index_data > start_index_data) ) {
            return 0.0;
          }

          size_t dims = dataset->getNrows();
          size_t sourceSize = source.getSize();
          double* ptrResult = result.getPointer();
          double time = 0.0;

          if (kernel_multTrans[0] == NULL) {
            OCLBasisType basis;
            basis.createMultTrans(dims, context, num_devices, device_ids, kernel_multTrans);
          }

          initOCLBuffers(level, index, mask, offset, dataset);

          cl_mem clSource[MAX_OCL_DEVICE_COUNT], clResult[MAX_OCL_DEVICE_COUNT];

          // create buffers for this execution
          for (size_t i = 0; i < num_devices; i++) {
            clSource[i] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * sourceSize, source.getPointer(), NULL);
            clResult[i] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(double) * result.getSize(), result.getPointer(), NULL);
          }

          // determine best fit
          size_t gpu_start_index_grid[MAX_OCL_DEVICE_COUNT];
          size_t gpu_end_index_grid[MAX_OCL_DEVICE_COUNT];

          for (size_t gpu_num = 0; gpu_num < num_devices; gpu_num++) {
            sg::parallel::PartitioningTool::getPartitionSegment(start_index_grid, end_index_grid, num_devices, gpu_num, &gpu_start_index_grid[gpu_num], &gpu_end_index_grid[gpu_num], OCL_SGPP_LOCAL_WORKGROUP_SIZE);
          }

          // set kernel arguments
          cl_uint clSourceSize = (cl_uint)sourceSize;
          cl_uint gpu_start_data = start_index_data;
          cl_uint gpu_end_data = end_index_data;

          for (size_t i = 0; i < num_devices; i++) {
            cl_uint gpu_start_grid = (cl_uint)gpu_start_index_grid[i];
            cl_uint gpu_end_grid = (cl_uint)gpu_end_index_grid[i];

            if (gpu_end_grid > gpu_start_grid) {
              if (clSetKernelArg(kernel_multTrans[i], 0, sizeof(cl_mem), &clLevel[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 1, sizeof(cl_mem), &clIndex[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 2, sizeof(cl_mem), &clMask[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 3, sizeof(cl_mem), &clOffset[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 4, sizeof(cl_mem), &clData[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 5, sizeof(cl_mem), &clSource[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 6, sizeof(cl_mem), &clResult[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 7, sizeof(cl_uint), &clSourceSize) || // sourceSize == number of entries in dataset
                  clSetKernelArg(kernel_multTrans[i], 8, sizeof(cl_uint), &gpu_start_data) ||
                  clSetKernelArg(kernel_multTrans[i], 9, sizeof(cl_uint), &gpu_end_data) != CL_SUCCESS ) {
                std::cout << "OCL Error: Failed to create kernel Args for kernel " << i << "!" << std::endl;
                return 0.0;
              }
            }
          }

          cl_event clTimings[MAX_OCL_DEVICE_COUNT];
          cl_event GPUDone[MAX_OCL_DEVICE_COUNT];

          // enqueue kernels
          size_t local = OCL_SGPP_LOCAL_WORKGROUP_SIZE;
          size_t active_devices = 0;

          for (size_t i = 0; i < num_devices; i++) {
            size_t rangeSize = gpu_end_index_grid[i] - gpu_start_index_grid[i];

            if (rangeSize > 0) {
              err = clEnqueueNDRangeKernel(command_queue[i], kernel_multTrans[i], 1, &gpu_start_index_grid[i], &rangeSize, &local, 0, NULL, &(clTimings[i]));

              if (active_devices != i) {
                std::cout << "OCL Error: Splitting up calculations for multiple GPUs is erroneous, only the last chunks may handle 0 entries. syncing will be not correct.";
              }

              active_devices++;

              if (err != CL_SUCCESS) {
                std::cout << "OCL Error: Failed to enqueue kernel command! Error Code: " << err << std::endl;
                return 0.0;
              }
            }
          }

          // read data back
          for (size_t i = 0; i < num_devices; i++) {
            size_t rangeSize = gpu_end_index_grid[i] - gpu_start_index_grid[i];

            if (rangeSize > 0) {
              size_t offset = gpu_start_index_grid[i];
              err = clEnqueueReadBuffer(command_queue[i], clResult[i], CL_FALSE, sizeof(double) * offset, sizeof(double) * rangeSize, &(ptrResult[offset]), 0, NULL, &(GPUDone[i]));

              if (err != CL_SUCCESS) {
                std::cout << "OCL Error: Failed to enqueue read buffer command (mult)! Error Code: " << err << std::endl;
                return 0.0;
              }
            }
          }

          // sync GPUs
          clWaitForEvents((cl_uint)active_devices, GPUDone);

          // determine kernel execution time
          for (size_t i = 0; i < num_devices; i++) {
            double tmpTime;
            cl_ulong startTime, endTime;
            startTime = endTime = 0;

            if (gpu_end_index_grid[i] > gpu_start_index_grid[i]) {
              err = clGetEventProfilingInfo(clTimings[i], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, NULL);

              if (err != CL_SUCCESS) {
                std::cout << "OCL Error: Failed to read start-time from command queue! Error Code: " << err << std::endl;
              }

              err = clGetEventProfilingInfo(clTimings[i], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, NULL);

              if (err != CL_SUCCESS) {
                std::cout << "OCL Error: Failed to read end-time from command queue! Error Code: " << err << std::endl;
              }
            }

            tmpTime = (double)(endTime - startTime);
            tmpTime *= 1e-9;

            if (tmpTime > time) {
              time = tmpTime;
            }
          }

          // clean up
          for (size_t i = 0; i < num_devices; i++) {
            clReleaseMemObject(clSource[i]);
            clReleaseMemObject(clResult[i]);

            if (gpu_end_index_grid[i] > gpu_start_index_grid[i]) {
              clReleaseEvent(clTimings[i]);
              clReleaseEvent(GPUDone[i]);
            }
          }

          return time;
        }

        void resetKernel() {
          // releaseKernelsAndPrograms();
          releaseGridBuffers();
        }

        static inline size_t getChunkGridPoints() {
          return 12;
        }

        static inline size_t getChunkDataPoints() {
          return OCL_SGPP_LOCAL_WORKGROUP_SIZE; /// TODO does this make sense?
        }

      private:
        void releaseGridBuffers() {
          for (size_t i = 0; i < num_devices; i++) {
            if (clLevel[i]) {
              clReleaseMemObject(clLevel[i]);
              clLevel[i] = NULL;
            }

            if (clIndex[i]) {
              clReleaseMemObject(clIndex[i]);
              clIndex[i] = NULL;
            }

            if (clMask[i]) {
              clReleaseMemObject(clMask[i]);
              clMask[i] = NULL;
            }

            if (clOffset[i]) {
              clReleaseMemObject(clOffset[i]);
              clOffset[i] = NULL;
            }
          }
        }

        void releaseDataBuffers() {
          for (size_t i = 0; i < num_devices; i++) {
            if (clData[i]) {
              clReleaseMemObject(clData[i]);
              clData[i] = NULL;
            }
          }
        }

        void releaseKernelsAndPrograms() {
          for (size_t i = 0; i < num_devices; i++) {
            if (kernel_mult[i]) {
              clReleaseKernel(kernel_mult[i]);
              kernel_mult[i] = NULL;
            }

            if (kernel_multTrans[i]) {
              clReleaseKernel(kernel_multTrans[i]);
              kernel_multTrans[i] = NULL;
            }
          }
        }
    };

  }
}

#endif // OCLKERNELIMPL_HPP

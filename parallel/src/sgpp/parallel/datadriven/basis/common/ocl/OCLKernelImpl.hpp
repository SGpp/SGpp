// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OCLKERNELIMPL_HPP
#define OCLKERNELIMPL_HPP

// include OpenCL
#include <CL/cl.h>

#include <string.h>

#include <sgpp/parallel/datadriven/basis/common/KernelBase.hpp>
#include <sgpp/parallel/datadriven/basis/common/ocl/OCLKernelImplBase.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    template<typename OCLBasisType>
    class OCLKernelImpl: public OCLKernelImplBase {
        double* pinnedGrid;
        double* pinnedTmp;
      public:
        OCLKernelImpl(): OCLKernelImplBase() {}

        inline void initOCLBuffers(
          SGPP::base::DataMatrix* level,
          SGPP::base::DataMatrix* index,
          SGPP::base::DataMatrix* mask,
          SGPP::base::DataMatrix* offset,
          SGPP::base::DataMatrix* dataset) {
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

        inline void initParams(SGPP::base::DataVector& grid,
                               SGPP::base::DataVector& tmp) {
          if (clPinnedGrid == NULL) {
            size_t mem_size = sizeof(double) * grid.getSize();
            clPinnedGrid = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, mem_size, NULL, NULL);

            for (size_t i = 0; i < num_devices; i++) {
              clDevGrid[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size, NULL, NULL);
            }

            pinnedGrid = (double*) clEnqueueMapBuffer(command_queue[0], clPinnedGrid, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE , 0, mem_size, 0, NULL, NULL, NULL);
          }

          for (size_t i = 0; i < grid.getSize(); i++) {
            pinnedGrid[i] = grid[i];
          }

          for (size_t i = 0; i < num_devices; i++) {
            clEnqueueWriteBuffer(command_queue[i], clDevGrid[i], CL_TRUE, 0, sizeof(double) * grid.getSize(), pinnedGrid, 0, NULL, NULL);
          }

          if (clPinnedTmp == NULL) {
            size_t mem_size = sizeof(double) * tmp.getSize();
            clPinnedTmp = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, mem_size, NULL, NULL);

            for (size_t i = 0; i < num_devices; i++) {
              clDevTmp[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size, NULL, NULL);
            }

            pinnedTmp = (double*) clEnqueueMapBuffer(command_queue[0], clPinnedTmp, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, mem_size, 0, NULL, NULL, NULL);
          }

          for (size_t i = 0; i < tmp.getSize(); i++) {
            pinnedTmp[i] = tmp[i];
          }

          for (size_t i = 0; i < num_devices; i++) {
            clEnqueueWriteBuffer(command_queue[i], clDevTmp[i], CL_TRUE, 0, sizeof(double) * tmp.getSize(), pinnedTmp, 0, NULL, NULL);
          }
        }

        double multImpl(
          SGPP::base::DataMatrix* level,
          SGPP::base::DataMatrix* index,
          SGPP::base::DataMatrix* mask,
          SGPP::base::DataMatrix* offset,
          SGPP::base::DataMatrix* dataset,
          SGPP::base::DataVector& alpha,
          SGPP::base::DataVector& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {
          // check if there is something to do at all
          if (!(end_index_grid > start_index_grid && end_index_data > start_index_data)) {
            return 0.0;
          }

          double time = 0.0;

          if (kernel_mult[0] == NULL) {
            OCLBasisType basis;
            size_t dims = dataset->getNrows();
            basis.createMult(dims, ocl_local_size, context, num_devices, device_ids, kernel_mult);
          }

          initOCLBuffers(level, index, mask, offset, dataset);
          initParams(alpha, result);

          // determine best fit
          size_t* gpu_start_index_data = new size_t[num_devices];
          size_t* gpu_end_index_data = new size_t[num_devices];

          for (size_t gpu_num = 0; gpu_num < num_devices; gpu_num++) {
            SGPP::parallel::PartitioningTool::getPartitionSegment(start_index_data, end_index_data, num_devices, gpu_num, &gpu_start_index_data[gpu_num], &gpu_end_index_data[gpu_num], ocl_local_size);
          }

          // set kernel arguments
          cl_uint clResultSize = (cl_uint)(result.getSize());
          cl_uint gpu_start_grid = (cl_uint)start_index_grid;
          cl_uint gpu_end_grid = (cl_uint)end_index_grid;

          for (size_t i = 0; i < num_devices; i++) {
            cl_uint gpu_start_data = (cl_uint)gpu_start_index_data[i];
            cl_uint gpu_end_data = (cl_uint)gpu_end_index_data[i];

            if (gpu_end_data > gpu_start_data) {
              if (clSetKernelArg(kernel_mult[i], 0, sizeof(cl_mem), &clLevel[i]) ||
                  clSetKernelArg(kernel_mult[i], 1, sizeof(cl_mem), &clIndex[i]) ||
                  clSetKernelArg(kernel_mult[i], 2, sizeof(cl_mem), &clMask[i]) || // only needed for masked version of modlinear
                  clSetKernelArg(kernel_mult[i], 3, sizeof(cl_mem), &clOffset[i]) || //only needed for masked version of modlinear
                  clSetKernelArg(kernel_mult[i], 4, sizeof(cl_mem), &clData[i]) ||
                  clSetKernelArg(kernel_mult[i], 5, sizeof(cl_mem), &clDevGrid[i]) ||
                  clSetKernelArg(kernel_mult[i], 6, sizeof(cl_mem), &clDevTmp[i]) ||
                  clSetKernelArg(kernel_mult[i], 7, sizeof(cl_uint), &clResultSize) || // resultsize == number of entries in dataset
                  clSetKernelArg(kernel_mult[i], 8, sizeof(cl_uint), &gpu_start_grid) ||
                  clSetKernelArg(kernel_mult[i], 9, sizeof(cl_uint), &gpu_end_grid) != CL_SUCCESS) {
                std::cout << "OCL Error: Failed to create kernel Args for mult!" << std::endl;
                return 0.0;
              }
            }
          }

          cl_event* clTimings = new cl_event[num_devices];
          cl_event* GPUDone = new cl_event[num_devices];

          // enqueue kernel
          size_t local = ocl_local_size;
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
              err = clEnqueueReadBuffer(command_queue[i], clDevTmp[i], CL_FALSE, sizeof(double) * offset, sizeof(double) * rangeSize, &(pinnedTmp[offset]), 0, NULL, &(GPUDone[i]));

              if (err != CL_SUCCESS) {
                std::cout << "OCL Error: Failed to enqueue read buffer command (mult)! Error Code: " << err << std::endl;
                return 0.0;
              }
            }
          }

          // sync GPUs
          clWaitForEvents((cl_uint)active_devices, GPUDone);

          for (size_t i = start_index_data; i < end_index_data; i++) {
            result[i] = pinnedTmp[i];
          }

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
            if (gpu_end_index_data[i] > gpu_start_index_data[i]) {
              clReleaseEvent(clTimings[i]);
              clReleaseEvent(GPUDone[i]);
            }
          }

          delete[] gpu_start_index_data;
          delete[] gpu_end_index_data;

          delete[] clTimings;
          delete[] GPUDone;

          return time;
        }

        double multTransposeImpl(
          SGPP::base::DataMatrix* level,
          SGPP::base::DataMatrix* index,
          SGPP::base::DataMatrix* mask,
          SGPP::base::DataMatrix* offset,
          SGPP::base::DataMatrix* dataset,
          SGPP::base::DataVector& source,
          SGPP::base::DataVector& result,
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
          double time = 0.0;

          if (kernel_multTrans[0] == NULL) {
            OCLBasisType basis;
            basis.createMultTrans(dims, ocl_local_size, context, num_devices, device_ids, kernel_multTrans);
          }

          initOCLBuffers(level, index, mask, offset, dataset);
          initParams(result, source);

          // determine best fit
          size_t* gpu_start_index_grid = new size_t[num_devices];
          size_t* gpu_end_index_grid = new size_t[num_devices];

          for (size_t gpu_num = 0; gpu_num < num_devices; gpu_num++) {
            SGPP::parallel::PartitioningTool::getPartitionSegment(start_index_grid, end_index_grid, num_devices, gpu_num, &gpu_start_index_grid[gpu_num], &gpu_end_index_grid[gpu_num], ocl_local_size);
          }

          // set kernel arguments
          cl_uint clSourceSize = (cl_uint)sourceSize;
          cl_uint gpu_start_data = (cl_uint)start_index_data;
          cl_uint gpu_end_data = (cl_uint)end_index_data;

          for (size_t i = 0; i < num_devices; i++) {
            cl_uint gpu_start_grid = (cl_uint)gpu_start_index_grid[i];
            cl_uint gpu_end_grid = (cl_uint)gpu_end_index_grid[i];

            if (gpu_end_grid > gpu_start_grid) {
              if (clSetKernelArg(kernel_multTrans[i], 0, sizeof(cl_mem), &clLevel[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 1, sizeof(cl_mem), &clIndex[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 2, sizeof(cl_mem), &clMask[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 3, sizeof(cl_mem), &clOffset[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 4, sizeof(cl_mem), &clData[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 5, sizeof(cl_mem), &clDevTmp[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 6, sizeof(cl_mem), &clDevGrid[i]) ||
                  clSetKernelArg(kernel_multTrans[i], 7, sizeof(cl_uint), &clSourceSize) || // sourceSize == number of entries in dataset
                  clSetKernelArg(kernel_multTrans[i], 8, sizeof(cl_uint), &gpu_start_data) ||
                  clSetKernelArg(kernel_multTrans[i], 9, sizeof(cl_uint), &gpu_end_data) != CL_SUCCESS ) {
                std::cout << "OCL Error: Failed to create kernel Args for kernel " << i << "!" << std::endl;
                return 0.0;
              }
            }
          }

          cl_event* clTimings = new cl_event[num_devices];
          cl_event* GPUDone = new cl_event[num_devices];

          // enqueue kernels
          size_t local = ocl_local_size;
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
              err = clEnqueueReadBuffer(command_queue[i], clDevGrid[i], CL_FALSE, sizeof(double) * offset, sizeof(double) * rangeSize, &(pinnedGrid[offset]), 0, NULL, &(GPUDone[i]));

              if (err != CL_SUCCESS) {
                std::cout << "OCL Error: Failed to enqueue read buffer command (multTrans)! Error Code: " << err << std::endl;
                return 0.0;
              }
            }
          }

          // sync GPUs
          clWaitForEvents((cl_uint)active_devices, GPUDone);

          for (size_t i = start_index_grid; i < end_index_grid; i++) {
            result[i] = pinnedGrid[i];
          }

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
            if (gpu_end_index_grid[i] > gpu_start_index_grid[i]) {
              clReleaseEvent(clTimings[i]);
              clReleaseEvent(GPUDone[i]);
            }
          }

          delete[] gpu_start_index_grid;
          delete[] gpu_end_index_grid;

          delete[] clTimings;
          delete[] GPUDone;

          return time;
        }
    };

  }
}

#endif // OCLKERNELIMPL_HPP
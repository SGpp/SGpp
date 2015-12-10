// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include "../OperationMultipleEvalStreamingModOCLMask/StreamingModOCLMaskKernelImpl.hpp"

#include <CL/cl.h>

#include <string.h>
#include <limits>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/LinearLoadBalancer.hpp>
#include <sgpp/base/opencl/OCLClonedBuffer.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/OCLStretchedBuffer.hpp>
#include "../OperationMultipleEvalStreamingModOCLMask/StreamingModOCLMaskKernelSourceBuilder.hpp"

#include <chrono>

namespace SGPP {
  namespace datadriven {

    template<typename real_type>
    class StreamingModOCLMaskKernelImpl {
      private:
        size_t dims;

        cl_int err;

        base::OCLClonedBuffer deviceData;
        base::OCLClonedBuffer deviceLevel;
        base::OCLClonedBuffer deviceIndex;
        base::OCLClonedBuffer deviceMask;
        base::OCLClonedBuffer deviceOffset;

        // use pinned memory (on host and device) to speed up data transfers from/to GPU
        base::OCLStretchedBuffer deviceGrid;
        real_type* hostGrid;
        base::OCLStretchedBuffer deviceTemp;
        real_type* hostTemp;

        cl_kernel* kernel_multTrans;
        cl_kernel* kernel_mult;

        double* deviceTimingsMult;
        double* deviceTimingsMultTranspose;

        StreamingModOCLMaskKernelSourceBuilder kernelSourceBuilder;
        std::shared_ptr<base::OCLManager> manager;
        std::shared_ptr<base::OCLOperationConfiguration> parameters;

        base::LinearLoadBalancer multLoadBalancer;
        base::LinearLoadBalancer multTransposeLoadBalancer;

      public:

        StreamingModOCLMaskKernelImpl(size_t dims, std::shared_ptr<base::OCLManager> manager,
                                      std::shared_ptr<base::OCLOperationConfiguration> parameters) :
          deviceData(manager), deviceLevel(manager), deviceIndex(manager), deviceMask(manager), deviceOffset(manager), deviceGrid(
            manager), deviceTemp(manager), kernelSourceBuilder(parameters, dims), manager(manager), parameters(
              parameters), multLoadBalancer(manager, this->parameters), multTransposeLoadBalancer(manager,
                  this->parameters) {

          this->dims = dims;
          this->deviceTimingsMult = new double[manager->num_devices];
          this->deviceTimingsMultTranspose = new double[manager->num_devices];

          for (size_t i = 0; i < manager->num_devices; i++) {
            //initialize with same timing to enforce equal problem sizes in the beginning
            this->deviceTimingsMult[i] = 1.0;
            this->deviceTimingsMultTranspose[i] = 1.0;
          }

          this->err = CL_SUCCESS;

          this->hostGrid = nullptr;
          this->hostTemp = nullptr;

          this->kernel_mult = new cl_kernel[manager->num_devices];
          this->kernel_multTrans = new cl_kernel[manager->num_devices];

          // initialize arrays
          for (size_t i = 0; i < manager->num_devices; i++) {
            this->kernel_mult[i] = nullptr;
            this->kernel_multTrans[i] = nullptr;
          }
        }

        ~StreamingModOCLMaskKernelImpl() {
          releaseDataBuffers();
          releaseGridBuffers();

          for (size_t i = 0; i < manager->num_devices; i++) {
            if (kernel_mult[i]) {
              clReleaseKernel(kernel_mult[i]);
              kernel_mult[i] = nullptr;
            }

            if (kernel_multTrans[i]) {
              clReleaseKernel(kernel_multTrans[i]);
              kernel_multTrans[i] = nullptr;
            }
          }

          //        // release command queue
          //        for (size_t i = 0; i < num_devices; i++) {
          //            if (command_queue[i]) {
          //                clReleaseCommandQueue(command_queue[i]);
          //            }
          //        }
          //
          //        // release context
          //        clReleaseContext(context);
          //
          //        delete[] this->command_queue;

          this->deviceData.freeBuffer();
          this->deviceLevel.freeBuffer();
          this->deviceIndex.freeBuffer();

          delete[] this->kernel_mult;
          delete[] this->kernel_multTrans;

          delete[] this->deviceTimingsMult;
          delete[] this->deviceTimingsMultTranspose;

        }

        void resetKernel() {
          //leads to a reallocation before next kernel execution
          releaseGridBuffers();
        }

        double mult(std::vector<real_type>& level, std::vector<real_type>& index, std::vector<real_type>& mask,
                    std::vector<real_type>& offset, size_t gridSize, std::vector<real_type>& dataset, size_t datasetSize,
                    std::vector<real_type>& alpha, std::vector<real_type>& result, const size_t start_index_grid,
                    const size_t end_index_grid, const size_t start_index_data, const size_t end_index_data) {

          // check if there is something to do at all
          if (!(end_index_grid > start_index_grid && end_index_data > start_index_data)) {
            return 0.0;
          }

          double time = 0.0;

          if (kernel_mult[0] == nullptr) {
            this->createMult(this->dims, (*parameters)["LOCAL_SIZE"].getUInt(), manager->context, manager->num_devices,
                             manager->device_ids, kernel_mult);
          }

          initOCLBuffers(level, index, mask, offset, gridSize, dataset, datasetSize);
          initParams(alpha, gridSize, result, datasetSize);

          // determine best fit
          size_t* gpu_start_index_data = new size_t[manager->num_devices];
          size_t* gpu_end_index_data = new size_t[manager->num_devices];

          //TODO: don't forget to set padding to DATA_BLOCKING * THREAD_BLOCK_SIZE

          multLoadBalancer.update(this->deviceTimingsMult);

          //        size_t dataBlockingSize = (*parameters)["KERNEL_DATA_BLOCKING_SIZE"].getUInt();
          size_t dataBlockingSize = 1;
          multLoadBalancer.getPartitionSegments(start_index_data, end_index_data, (*parameters)["LOCAL_SIZE"].getUInt(),
                                                gpu_start_index_data, gpu_end_index_data);

          // set kernel arguments
          cl_uint clResultSize = (cl_uint) datasetSize;
          cl_uint gpu_start_grid = (cl_uint) start_index_grid;
          cl_uint gpu_end_grid = (cl_uint) end_index_grid;
          //    std::cout << "start grid: " << gpu_start_grid << " end grid: " << gpu_end_grid << std::endl;

          for (size_t i = 0; i < manager->num_devices; i++) {
            cl_uint gpu_start_data = static_cast<cl_uint>(gpu_start_index_data[i])
                                     / static_cast<cl_uint>(dataBlockingSize);
            cl_uint gpu_end_data = static_cast<cl_uint>(gpu_end_index_data[i]) / static_cast<cl_uint>(dataBlockingSize);
            //      std::cout << "device: " << i << " start data: " << gpu_start_data << " end data: " << gpu_end_data << std::endl;

            if (gpu_end_data > gpu_start_data) {
              err = clSetKernelArg(kernel_mult[i], 0, sizeof(cl_mem), this->deviceLevel.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_mult[i], 1, sizeof(cl_mem), this->deviceIndex.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_mult[i], 2, sizeof(cl_mem), this->deviceMask.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_mult[i], 3, sizeof(cl_mem), this->deviceOffset.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_mult[i], 4, sizeof(cl_mem), this->deviceData.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_mult[i], 5, sizeof(cl_mem), this->deviceGrid.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_mult[i], 6, sizeof(cl_mem), this->deviceTemp.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_mult[i], 7, sizeof(cl_uint), &clResultSize); // resultsize == number of entries in dataset

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_mult[i], 8, sizeof(cl_uint), &gpu_start_grid);

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_mult[i], 9, sizeof(cl_uint), &gpu_end_grid);

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }
            }
          }

          cl_event* clTimings = new cl_event[manager->num_devices];

          // enqueue kernel
          size_t local = (*parameters)["LOCAL_SIZE"].getUInt();
          size_t active_devices = 0;

          std::chrono::time_point<std::chrono::system_clock> start, end;
          start = std::chrono::system_clock::now();

          for (size_t i = 0; i < manager->num_devices; i++) {
            size_t rangeSize = (gpu_end_index_data[i] / dataBlockingSize)
                               - (gpu_start_index_data[i] / dataBlockingSize);

            if (rangeSize > 0) {
              err = clEnqueueNDRangeKernel(manager->command_queue[i], kernel_mult[i], 1, &gpu_start_index_data[i],
                                           &rangeSize, &local, 0, nullptr, &(clTimings[i]));

              if (active_devices != i) {
                std::stringstream errorString;
                errorString
                    << "OCL Error: Multiple GPUs is erroneous, because only the last chunks may handle 0 entries. Synchronization will be not correct."
                    << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              active_devices++;

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to enqueue kernel command! Error code: " << err << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }
            }
          }

          deviceTemp.readFromBuffer(gpu_start_index_data, gpu_end_index_data);

          end = std::chrono::system_clock::now();
          //        std::chrono::duration<double> elapsed_seconds = end - start;

          for (size_t i = start_index_data; i < end_index_data; i++) {
            result[i] = hostTemp[i];
          }

          // determine kernel execution time
          for (size_t i = 0; i < manager->num_devices; i++) {
            double tmpTime;
            cl_ulong startTime, endTime;
            startTime = endTime = 0;

            if (gpu_end_index_data[i] > gpu_start_index_data[i]) {
              err = clGetEventProfilingInfo(clTimings[i],
                                            CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, nullptr);

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString
                    << "OCL Error: Failed to read start-time from command queue (or crash in mult)! Error code: "
                    << err << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clGetEventProfilingInfo(clTimings[i],
                                            CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, nullptr);

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to read end-time from command queue! Error code: " << err
                            << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }
            }

            tmpTime = (double) (endTime - startTime);
            tmpTime *= 1e-9;
            this->deviceTimingsMult[i] = tmpTime;

            if (tmpTime > time) {
              time = tmpTime;
            }
          }

          // clean up
          for (size_t i = 0; i < manager->num_devices; i++) {
            if (gpu_end_index_data[i] > gpu_start_index_data[i]) {
              clReleaseEvent(clTimings[i]);
            }
          }

          delete[] gpu_start_index_data;
          delete[] gpu_end_index_data;

          delete[] clTimings;

          return time;
        }

        double multTranspose(std::vector<real_type>& level, std::vector<real_type>& index, std::vector<real_type>& mask,
                             std::vector<real_type>& offset, size_t gridSize, std::vector<real_type>& dataset, size_t datasetSize,
                             std::vector<real_type>& source, std::vector<real_type>& result, const size_t start_index_grid,
                             const size_t end_index_grid, const size_t start_index_data, const size_t end_index_data) {

          // check if there is something to do at all
          if (!(end_index_grid > start_index_grid && end_index_data > start_index_data)) {
            return 0.0;
          }

          double time = 0.0;

          if (kernel_multTrans[0] == nullptr) {
            this->createMultTrans(this->dims, (*parameters)["LOCAL_SIZE"].getUInt(), manager->context, manager->num_devices,
                                  manager->device_ids, kernel_multTrans);
          }

          initOCLBuffers(level, index, mask, offset, gridSize, dataset, datasetSize);
          initParams(result, gridSize, source, datasetSize);

          // determine best fit
          size_t* gpu_start_index_grid = new size_t[manager->num_devices];
          size_t* gpu_end_index_grid = new size_t[manager->num_devices];

          multTransposeLoadBalancer.update(this->deviceTimingsMultTranspose);
          multTransposeLoadBalancer.getPartitionSegments(start_index_grid, end_index_grid,
              (*parameters)["LOCAL_SIZE"].getUInt(), gpu_start_index_grid, gpu_end_index_grid);

          // set kernel arguments
          cl_uint clSourceSize = (cl_uint) datasetSize;
          cl_uint gpu_start_data = (cl_uint) start_index_data;
          cl_uint gpu_end_data = (cl_uint) end_index_data;

          //    std::cout << "start data: " << gpu_start_data << " end data: " << gpu_end_data << std::endl;

          for (size_t i = 0; i < manager->num_devices; i++) {
            cl_uint gpu_start_grid = (cl_uint) gpu_start_index_grid[i];
            cl_uint gpu_end_grid = (cl_uint) gpu_end_index_grid[i];
            //      std::cout << "device: " << i << " start grid: " << gpu_start_grid << " end grid: " << gpu_end_grid << std::endl;

            if (gpu_end_grid > gpu_start_grid) {
              err = clSetKernelArg(kernel_multTrans[i], 0, sizeof(cl_mem), this->deviceLevel.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_multTrans[i], 1, sizeof(cl_mem), this->deviceIndex.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_multTrans[i], 2, sizeof(cl_mem), this->deviceMask.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_multTrans[i], 3, sizeof(cl_mem), this->deviceOffset.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_multTrans[i], 4, sizeof(cl_mem), this->deviceData.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_multTrans[i], 5, sizeof(cl_mem), this->deviceTemp.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_multTrans[i], 6, sizeof(cl_mem), this->deviceGrid.getBuffer(i));

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_multTrans[i], 7, sizeof(cl_uint), &clSourceSize); // sourceSize == number of entries in dataset

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_multTrans[i], 8, sizeof(cl_uint), &gpu_start_data);

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clSetKernelArg(kernel_multTrans[i], 9, sizeof(cl_uint), &gpu_end_data);

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }
            }
          }

          cl_event* clTimings = new cl_event[manager->num_devices];

          // enqueue kernels
          size_t local = (*parameters)["LOCAL_SIZE"].getUInt();
          size_t active_devices = 0;

          for (size_t i = 0; i < manager->num_devices; i++) {
            size_t rangeSize = gpu_end_index_grid[i] - gpu_start_index_grid[i];

            if (rangeSize > 0) {
              //      std::cout << "enqueuing device: " << i << std::endl;
              err = clEnqueueNDRangeKernel(manager->command_queue[i], kernel_multTrans[i], 1,
                                           &gpu_start_index_grid[i], &rangeSize, &local, 0, nullptr, &(clTimings[i]));

              if (active_devices != i) {
                std::stringstream errorString;
                errorString
                    << "OCL Error: Multiple GPUs is erroneous, because only the last chunks may handle 0 entries. Synchronization will be not correct."
                    << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              active_devices++;

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to enqueue kernel command! Error code: " << err << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }
            }
          }

          deviceGrid.readFromBuffer(gpu_start_index_grid, gpu_end_index_grid);

          for (size_t i = start_index_grid; i < end_index_grid; i++) {
            result[i] = hostGrid[i];
          }

          // determine kernel execution time
          for (size_t i = 0; i < manager->num_devices; i++) {
            double tmpTime;
            cl_ulong startTime, endTime;
            startTime = endTime = 0;

            if (gpu_end_index_grid[i] > gpu_start_index_grid[i]) {
              err = clGetEventProfilingInfo(clTimings[i],
                                            CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, nullptr);

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString
                    << "OCL Error: Failed to read start-time from command queue (or crash in multTranspose)! Error code: "
                    << err << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }

              err = clGetEventProfilingInfo(clTimings[i],
                                            CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, nullptr);

              if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to read end-time from command queue! Error code: " << err
                            << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
              }
            }

            tmpTime = (double) (endTime - startTime);
            tmpTime *= 1e-9;
            this->deviceTimingsMultTranspose[i] = tmpTime;

            if (tmpTime > time) {
              time = tmpTime;
            }
          }

          // clean up
          for (size_t i = 0; i < manager->num_devices; i++) {
            if (gpu_end_index_grid[i] > gpu_start_index_grid[i]) {
              clReleaseEvent(clTimings[i]);
            }
          }

          delete[] gpu_start_index_grid;
          delete[] gpu_end_index_grid;

          delete[] clTimings;

          return time;
        }
      private:

        void createMultTrans(size_t dims, size_t local_workgroup_size, cl_context context, size_t num_devices,
                             cl_device_id* device_ids, cl_kernel* kernel) {
          std::string program_src = kernelSourceBuilder.generateSourceMultTrans();
          manager->buildKernel(program_src, "multTransOCLMask", context, num_devices, device_ids, kernel);
        }

        void createMult(size_t dims, size_t local_workgroup_size, cl_context context, size_t num_devices,
                        cl_device_id* device_ids, cl_kernel* kernel) {
          std::string program_src = kernelSourceBuilder.generateSourceMult();
          manager->buildKernel(program_src, "multOCLMask", context, num_devices, device_ids, kernel);
        }

        void releaseGridBuffers() {
          this->deviceLevel.freeBuffer();
          this->deviceIndex.freeBuffer();
          this->deviceMask.freeBuffer();
          this->deviceOffset.freeBuffer();
          this->deviceGrid.freeBuffer();
        }

        void releaseDataBuffers() {
          this->deviceData.freeBuffer();
          this->deviceTemp.freeBuffer();
        }

        void initOCLBuffers(std::vector<real_type>& level, std::vector<real_type>& index, std::vector<real_type>& mask,
                            std::vector<real_type>& offset, size_t gridSize, std::vector<real_type>& dataset, size_t datasetSize) {
          if (!this->deviceLevel.isInitialized()) {
            this->deviceLevel.initializeBuffer(level.data(), sizeof(real_type), gridSize * this->dims);
          }

          if (!this->deviceIndex.isInitialized()) {
            this->deviceIndex.initializeBuffer(index.data(), sizeof(real_type), gridSize * this->dims);
          }

          if (!this->deviceMask.isInitialized()) {
            this->deviceMask.initializeBuffer(mask.data(), sizeof(real_type), gridSize * this->dims);
          }

          if (!this->deviceOffset.isInitialized()) {
            this->deviceOffset.initializeBuffer(offset.data(), sizeof(real_type), gridSize * this->dims);
          }

          if (!this->deviceData.isInitialized()) {
            this->deviceData.initializeBuffer(dataset.data(), sizeof(real_type), datasetSize * this->dims);
          }
        }

        void initParams(std::vector<real_type>& grid, size_t gridSize, std::vector<real_type>& tmp, size_t datasetSize) {
          if (!this->deviceGrid.isInitialized()) {
            this->deviceGrid.initializeBuffer(sizeof(real_type), gridSize * this->dims);
            this->hostGrid = (real_type*) this->deviceGrid.getMappedHostBuffer();
          }

          for (size_t i = 0; i < gridSize; i++) {
            this->hostGrid[i] = grid[i];
          }

          deviceGrid.writeToBuffer();

          if (!this->deviceTemp.isInitialized()) {
            this->deviceTemp.initializeBuffer(sizeof(real_type), datasetSize * this->dims);
            this->hostTemp = (real_type*) this->deviceTemp.getMappedHostBuffer();
          }

          for (size_t i = 0; i < datasetSize; i++) {
            this->hostTemp[i] = tmp[i];
          }

          deviceTemp.writeToBuffer();
        }

    };

  }
}

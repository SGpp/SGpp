// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <CL/cl.h>
#include <omp.h>

#include <string>
#include <limits>
#include <chrono>
#include <vector>

#include "sgpp/base/opencl/OCLBufferWrapperSD.hpp"
#include "sgpp/globaldef.hpp"
#include "sgpp/base/opencl/LinearLoadBalancerMultiPlatform.hpp"
#include "sgpp/base/opencl/OCLManagerMultiPlatform.hpp"
#include "sgpp/base/opencl/OCLStretchedBuffer.hpp"
#include "SourceBuilderMult.hpp"

namespace SGPP {
namespace datadriven {
namespace StreamingModOCLMaskMultiPlatform {

template <typename real_type>
class KernelMult {
 private:
  std::shared_ptr<base::OCLDevice> device;

  size_t dims;

  cl_int err;

  base::OCLBufferWrapperSD<real_type> deviceLevel;
  base::OCLBufferWrapperSD<real_type> deviceIndex;
  base::OCLBufferWrapperSD<real_type> deviceMask;
  base::OCLBufferWrapperSD<real_type> deviceOffset;
  base::OCLBufferWrapperSD<real_type> deviceAlpha;

  base::OCLBufferWrapperSD<real_type> deviceData;

  base::OCLBufferWrapperSD<real_type> deviceResultData;

  cl_kernel kernelMult;

  double deviceTimingMult;

  SourceBuilderMult<real_type> kernelSourceBuilder;
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  //    std::shared_ptr<base::OCLOperationConfiguration> parameters;
  json::Node &kernelConfiguration;

  std::shared_ptr<base::QueueLoadBalancer> queueLoadBalancerMult;

  bool verbose;

  size_t localSize;
  size_t dataBlockingSize;
  size_t scheduleSize;
  size_t totalBlockSize;

 public:
  KernelMult(std::shared_ptr<base::OCLDevice> device, size_t dims,
             std::shared_ptr<base::OCLManagerMultiPlatform> manager,
             json::Node &kernelConfiguration,
             std::shared_ptr<base::QueueLoadBalancer> queueBalancerMult)
      : device(device),
        dims(dims),
        err(CL_SUCCESS),
        deviceLevel(device),
        deviceIndex(device),
        deviceMask(device),
        deviceOffset(device),
        deviceAlpha(device),
        deviceData(device),
        deviceResultData(device),
        kernelMult(nullptr),
        kernelSourceBuilder(device, kernelConfiguration, dims),
        manager(manager),
        kernelConfiguration(kernelConfiguration),
        queueLoadBalancerMult(queueBalancerMult) {
    // initialize with same timing to enforce equal problem sizes in the
    // beginning
    this->deviceTimingMult = 1.0;
    this->verbose = kernelConfiguration["VERBOSE"].getBool();

    localSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    dataBlockingSize = kernelConfiguration["KERNEL_DATA_BLOCK_SIZE"].getUInt();
    scheduleSize = kernelConfiguration["KERNEL_SCHEDULE_SIZE"].getUInt();
    totalBlockSize = localSize * dataBlockingSize;
  }

  ~KernelMult() {
    if (this->kernelMult != nullptr) {
      clReleaseKernel(this->kernelMult);
      this->kernelMult = nullptr;
    }
  }

  void resetKernel() {
    //    releaseGridBuffers();
    //    releaseDataBuffers();
    //    releaseDatasetResultBuffer();
  }

  double mult(std::vector<real_type> &level, std::vector<real_type> &index,
              std::vector<real_type> &mask, std::vector<real_type> &offset,
              std::vector<real_type> &dataset, std::vector<real_type> &alpha,
              std::vector<real_type> &result, const size_t start_index_grid,
              const size_t end_index_grid, const size_t start_index_data,
              const size_t end_index_data) {
    // check if there is something to do at all
    if (!(end_index_grid > start_index_grid && end_index_data > start_index_data)) {
      return 0.0;
    }

    if (this->kernelMult == nullptr) {
      std::string program_src = kernelSourceBuilder.generateSource();
      this->kernelMult =
          manager->buildKernel(program_src, device, kernelConfiguration, "multOCLMask");
    }

    if (!deviceLevel.isInitialized()) {
      initGridBuffers(level, index, mask, offset, alpha, start_index_grid, end_index_grid);
    }

    this->deviceTimingMult = 0.0;

    while (true) {
      size_t kernelStartData;
      size_t kernelEndData;

      // set kernel arguments
      size_t kernelStartGrid = start_index_grid;
      size_t kernelEndGrid = end_index_grid;

      bool segmentAvailable = queueLoadBalancerMult->getNextSegment(scheduleSize, totalBlockSize,
                                                                    kernelStartData, kernelEndData);
      if (!segmentAvailable) {
        break;
      }

      size_t rangeSizeUnblocked = kernelEndData - kernelStartData;

      if (verbose) {
        std::cout << "device: " << device->deviceId << " kernel from: " << kernelStartData
                  << " to: " << kernelEndData << " -> range: " << rangeSizeUnblocked << std::endl;
      }

      initDatasetBuffers(dataset, kernelStartData, kernelEndData);
      initDatasetResultBuffers(kernelStartData, kernelEndData);

      clFinish(device->commandQueue);
      //            std::cout << "wrote to device: " << device->deviceId << ""
      //            << std::endl;

      size_t rangeSizeBlocked =
          (kernelEndData / dataBlockingSize) - (kernelStartData / dataBlockingSize);

      if (rangeSizeBlocked > 0) {
        err = clSetKernelArg(this->kernelMult, 0, sizeof(cl_mem), this->deviceLevel.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 1, sizeof(cl_mem), this->deviceIndex.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 2, sizeof(cl_mem), this->deviceMask.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 3, sizeof(cl_mem), this->deviceOffset.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 4, sizeof(cl_mem), this->deviceData.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 5, sizeof(cl_mem), this->deviceAlpha.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }
        err =
            clSetKernelArg(this->kernelMult, 6, sizeof(cl_mem), this->deviceResultData.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 7, sizeof(cl_uint),
                             &rangeSizeUnblocked);  // resultsize == number of entries in dataset
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 8, sizeof(cl_uint), &kernelStartGrid);
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 9, sizeof(cl_uint), &kernelEndGrid);
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }

        cl_event clTiming = nullptr;

        // enqueue kernel
        //                size_t localSize =
        //                kernelConfiguration["LOCAL_SIZE"].getUInt();

        //                std::cout << "commandQueue: " << device->commandQueue
        //                << std::endl;

        char deviceName[128] = {0};
        cl_uint err;
        err = clGetDeviceInfo(device->deviceId, CL_DEVICE_NAME, 128 * sizeof(char), &deviceName,
                              nullptr);

        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to read the device name for device: "
                      << device->deviceId << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }

        //                std::cout << "OCL Info: detected device, name: \"" <<
        //                deviceName << "\"" << std::endl;

        err = clEnqueueNDRangeKernel(device->commandQueue, this->kernelMult, 1, 0,
                                     &rangeSizeBlocked, &localSize, 0, nullptr, &clTiming);

        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to enqueue kernel command! Error code: " << err
                      << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }

        clFinish(device->commandQueue);
        //                std::cout << "executed kernel: " << device->deviceId
        //                << "" << std::endl;

        deviceResultData.readFromBuffer();

        //                std::cout << "read from device: " << device->deviceId
        //                << "" << std::endl;
        clFinish(device->commandQueue);

        std::vector<real_type> &hostTemp = deviceResultData.getHostPointer();
        size_t deviceIndex = 0;
        for (size_t i = 0; i < rangeSizeUnblocked; i++) {
          //                    std::cout << "resultDevice[" << deviceIndex <<
          //                    "] = " << hostTemp[deviceIndex] << std::endl;
          result[kernelStartData + i] = hostTemp[deviceIndex];
          deviceIndex += 1;
        }

        // determine kernel execution time
        cl_ulong startTime = 0;
        cl_ulong endTime = 0;

        err = clGetEventProfilingInfo(clTiming, CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
                                      &startTime, nullptr);

        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to read start-time from command "
                         "queue (or crash in mult)! Error code: " << err << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }

        err = clGetEventProfilingInfo(clTiming, CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
                                      &endTime, nullptr);

        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to read end-time from command "
                         "queue! Error code: " << err << std::endl;
          throw SGPP::base::operation_exception(errorString.str());
        }

        clReleaseEvent(clTiming);

        double time = 0.0;
        time = static_cast<double>(endTime - startTime);
        time *= 1e-9;

        if (verbose) {
          std::cout << "device: " << device->deviceId << " duration: " << time << std::endl;
        }

        this->deviceTimingMult += time;
      }
    }

    return this->deviceTimingMult;
  }

 private:
  //  void releaseGridBuffers() {
  //    this->deviceLevel.freeBuffer();
  //    this->deviceIndex.freeBuffer();
  //    this->deviceMask.freeBuffer();
  //    this->deviceOffset.freeBuffer();
  //    this->deviceAlpha.freeBuffer();
  //  }
  //
  //  void releaseDataBuffers() {
  //    this->deviceData.freeBuffer();
  //    this->deviceResultData.freeBuffer();
  //  }
  //
  //  void releaseDatasetResultBuffer() { this->deviceResultData.freeBuffer(); }

  void initGridBuffers(std::vector<real_type> &level, std::vector<real_type> &index,
                       std::vector<real_type> &mask, std::vector<real_type> &offset,
                       std::vector<real_type> &alpha, size_t kernelStartGrid,
                       size_t kernelEndGrid) {
    deviceLevel.intializeTo(level, dims, kernelStartGrid, kernelEndGrid);
    deviceIndex.intializeTo(index, dims, kernelStartGrid, kernelEndGrid);
    deviceMask.intializeTo(mask, dims, kernelStartGrid, kernelEndGrid);
    deviceOffset.intializeTo(offset, dims, kernelStartGrid, kernelEndGrid);
    deviceAlpha.intializeTo(alpha, 1, kernelStartGrid, kernelEndGrid);
  }

  void initDatasetBuffers(std::vector<real_type> &dataset, size_t kernelStartData,
                          size_t kernelEndData) {
    deviceData.intializeTo(dataset, dims, kernelStartData, kernelEndData, true);
  }

  void initDatasetResultBuffers(size_t kernelStartData, size_t kernelEndData) {
    size_t range = kernelEndData - kernelStartData;

    std::vector<real_type> zeros(range);
    for (size_t i = 0; i < range; i++) {
      zeros[i] = 0.0;
    }

    deviceResultData.intializeTo(zeros, 1, 0, range);
  }
};
}  // namespace StreamingModOCLMaskMultiPlatform
}  // namespace datadriven
}  // namespace SGPP

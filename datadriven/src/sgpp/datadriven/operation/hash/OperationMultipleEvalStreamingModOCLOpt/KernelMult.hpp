// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <CL/cl.h>

#include <limits>
#include <string>
#include <vector>

#include <sgpp/base/tools/QueueLoadBalancerOpenMP.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLOpt/SourceBuilderMult.hpp>
#include <sgpp/base/opencl/OCLBufferWrapperSD.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {
namespace StreamingModOCLOpt {

template <typename T> class KernelMult {
private:
  std::shared_ptr<base::OCLDevice> device;

  size_t dims;

  cl_int err;

  base::OCLBufferWrapperSD<T> deviceData;
  base::OCLBufferWrapperSD<T> deviceLevel;
  base::OCLBufferWrapperSD<T> deviceIndex;

  base::OCLBufferWrapperSD<T> deviceAlpha;
  base::OCLBufferWrapperSD<T> deviceResultData;

  cl_kernel kernelMult;

  SourceBuilderMult<T> kernelSourceBuilder;

  std::shared_ptr<base::OCLManagerMultiPlatform> manager;

  double deviceTimingMult;

  json::Node &kernelConfiguration;

  std::shared_ptr<base::QueueLoadBalancerOpenMP> queueLoadBalancerMult;

  bool verbose;

  size_t localSize;
  size_t dataBlockSize;
  size_t scheduleSize;
  size_t totalBlockSize;

public:
  KernelMult(std::shared_ptr<base::OCLDevice> device, size_t dims,
             std::shared_ptr<base::OCLManagerMultiPlatform> manager,
             json::Node &kernelConfiguration,
             std::shared_ptr<base::QueueLoadBalancerOpenMP> queueBalancerMult)
      : device(device), dims(dims), err(CL_SUCCESS), deviceData(device),
        deviceLevel(device), deviceIndex(device), deviceAlpha(device),
        deviceResultData(device), kernelMult(nullptr),
        kernelSourceBuilder(device, kernelConfiguration, dims),
        manager(manager), deviceTimingMult(0.0),
        kernelConfiguration(kernelConfiguration),
        queueLoadBalancerMult(queueBalancerMult) {
    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") ==
            0 &&
        dims > kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt()) {
      std::stringstream errorString;
      errorString << "OCL Error: setting \"KERNEL_DATA_STORE\" to \"register\" "
                     "requires value of "
                     "\"KERNEL_MAX_DIM_UNROLL\" to be greater than the "
                     "dimension of the data "
                     "set, was set to "
                  << kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt()
                  << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    this->verbose = kernelConfiguration["VERBOSE"].getBool();

    localSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    dataBlockSize = kernelConfiguration["KERNEL_DATA_BLOCK_SIZE"].getUInt();
    scheduleSize = kernelConfiguration["KERNEL_SCHEDULE_SIZE"].getUInt();
    totalBlockSize = dataBlockSize * localSize;
  }

  ~KernelMult() {
    if (kernelMult != nullptr) {
      clReleaseKernel(kernelMult);
      this->kernelMult = nullptr;
    }
  }

  void resetKernel() {
    // leads to a reallocation before next kernel execution
    //    releaseGridBuffers();
  }

  double mult(std::vector<T> &level, std::vector<T> &index,
              std::vector<T> &dataset, std::vector<T> &alpha,
              std::vector<T> &result, const size_t start_index_grid,
              const size_t end_index_grid, const size_t start_index_data,
              const size_t end_index_data) {
    if (verbose) {
#pragma omp critical(StreamingOCLMultiPlatformKernelMultTranspose)
      {
        std::cout << "entering mult, device: " << device->deviceName << " ("
                  << device->deviceId << ")" << std::endl;
      }
    }

    // check if there is something to do at all
    if (!(end_index_grid > start_index_grid &&
          end_index_data > start_index_data)) {
      return 0.0;
    }

    if (this->kernelMult == nullptr) {
      std::string program_src = kernelSourceBuilder.generateSource();
      this->kernelMult = manager->buildKernel(program_src, device,
                                              kernelConfiguration, "multOCL");
    }

    initGridBuffers(level, index, alpha, start_index_grid, end_index_grid);

    this->deviceTimingMult = 0.0;

// for slow devices to catch up
#pragma omp barrier

    while (true) {
      size_t kernelStartData = 0;
      size_t kernelEndData = 0;

      // set kernel arguments
      size_t kernelStartGrid = start_index_grid;
      size_t kernelEndGrid = end_index_grid;

      bool segmentAvailable = queueLoadBalancerMult->getNextSegment(
          scheduleSize, kernelStartData, kernelEndData);
      if (!segmentAvailable) {
        break;
      }

      size_t rangeSize = kernelEndData - kernelStartData;
      size_t rangeSizeAfterBlocking =
          (kernelEndData / dataBlockSize) - (kernelStartData / dataBlockSize);

      if (verbose) {
#pragma omp critical(StreamingOCLMultiPlatformKernelMult)
        {
          std::cout << "device: " << device->deviceName << " ("
                    << device->deviceId << ") "
                    << " kernel from: " << kernelStartData
                    << " to: " << kernelEndData << " -> range: " << rangeSize
                    << " (with blocking: " << rangeSizeAfterBlocking << ")"
                    << std::endl;
        }
      }

      initDatasetBuffers(dataset, kernelStartData, kernelEndData);
      initDatasetResultBuffers(kernelStartData, kernelEndData);

      clFinish(device->commandQueue);
      //            std::cout << "wrote to device: " << device->deviceId << ""
      //            << std::endl;

      if (rangeSize > 0) {
        err = clSetKernelArg(this->kernelMult, 0, sizeof(cl_mem),
                             this->deviceLevel.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 1, sizeof(cl_mem),
                             this->deviceIndex.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 2, sizeof(cl_mem),
                             this->deviceData.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 3, sizeof(cl_mem),
                             this->deviceAlpha.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 4, sizeof(cl_mem),
                             this->deviceResultData.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(
            this->kernelMult, 5, sizeof(cl_int),
            &rangeSize); // resultsize == number of entries in dataset
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(this->kernelMult, 6, sizeof(cl_int),
                             &kernelStartGrid);
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err =
            clSetKernelArg(this->kernelMult, 7, sizeof(cl_int), &kernelEndGrid);
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }

        cl_event clTiming = nullptr;

        // enqueue kernel

        //                std::cout << "commandQueue: " << device->commandQueue
        //                << std::endl;

        char deviceName[128] = {0};
        cl_uint err;
        err = clGetDeviceInfo(device->deviceId, CL_DEVICE_NAME,
                              128 * sizeof(char), &deviceName, nullptr);

        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to read the device name for device: "
              << device->deviceId << std::endl;
          throw base::operation_exception(errorString.str());
        }

        //                std::cout << "OCL Info: detected device, name: \"" <<
        //                deviceName << "\"" << std::endl;

        err = clEnqueueNDRangeKernel(device->commandQueue, this->kernelMult, 1,
                                     0, &rangeSizeAfterBlocking, &localSize, 0,
                                     nullptr, &clTiming);

        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to enqueue kernel command! Error code: "
              << err << std::endl;
          throw base::operation_exception(errorString.str());
        }

        clFinish(device->commandQueue);
        //                std::cout << "executed kernel: " << device->deviceId
        //                << "" << std::endl;

        deviceResultData.readFromBuffer();

        //                std::cout << "read from device: " << device->deviceId
        //                << "" << std::endl;
        clFinish(device->commandQueue);

        std::vector<T> &hostTemp = deviceResultData.getHostPointer();
        size_t deviceIndex = 0;
        for (size_t i = 0; i < rangeSize; i++) {
          result[kernelStartData + i] = hostTemp[deviceIndex];
          deviceIndex += 1;
        }

        // determine kernel execution time
        cl_ulong startTime = 0;
        cl_ulong endTime = 0;

        err = clGetEventProfilingInfo(clTiming, CL_PROFILING_COMMAND_START,
                                      sizeof(cl_ulong), &startTime, nullptr);

        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to read start-time from command "
                         "queue (or crash in mult)! Error code: "
                      << err << std::endl;
          throw base::operation_exception(errorString.str());
        }

        err = clGetEventProfilingInfo(clTiming, CL_PROFILING_COMMAND_END,
                                      sizeof(cl_ulong), &endTime, nullptr);

        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to read end-time from command "
                         "queue! Error code: "
                      << err << std::endl;
          throw base::operation_exception(errorString.str());
        }

        clReleaseEvent(clTiming);

        double time = 0.0;
        time = static_cast<double>(endTime - startTime);
        time *= 1e-9;

        if (verbose) {
#pragma omp critical(StreamingOCLMultiPlatformKernelMult)
          {
            std::cout << "device: " << device->deviceName << " ("
                      << device->deviceId << ") "
                      << "duration: " << time << std::endl;
          }
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
  //    this->deviceAlpha.freeBuffer();
  //  }
  //
  //  void releaseDataBuffers() {
  //    this->deviceData.freeBuffer();
  //    this->deviceResultData.freeBuffer();
  //  }
  //
  //  void releaseDatasetResultBuffer() { this->deviceResultData.freeBuffer(); }

  void initGridBuffers(std::vector<T> &level, std::vector<T> &index,
                       std::vector<T> &alpha, size_t kernelStartGrid,
                       size_t kernelEndGrid) {
    deviceLevel.intializeTo(level, dims, kernelStartGrid, kernelEndGrid);
    deviceIndex.intializeTo(index, dims, kernelStartGrid, kernelEndGrid);
    deviceAlpha.intializeTo(alpha, 1, kernelStartGrid, kernelEndGrid);
  }

  void initDatasetBuffers(std::vector<T> &dataset, size_t kernelStartData,
                          size_t kernelEndData) {
    deviceData.intializeTo(dataset, dims, kernelStartData, kernelEndData, true);
  }

  void initDatasetResultBuffers(size_t kernelStartData, size_t kernelEndData) {
    size_t range = kernelEndData - kernelStartData;

    // TODO(pfandedd): use result parameter in mult
    std::vector<T> zeros(range);
    for (size_t i = 0; i < range; i++) {
      zeros[i] = 0.0;
    }

    deviceResultData.intializeTo(zeros, 1, 0, range);
  }
};
} // namespace StreamingModOCLOpt
} // namespace datadriven
} // namespace sgpp

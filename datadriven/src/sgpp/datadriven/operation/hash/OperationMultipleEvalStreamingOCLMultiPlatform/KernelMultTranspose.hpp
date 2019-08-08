// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <CL/cl.h>

#include <chrono>
#include <limits>
#include <string>
#include <vector>

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/SourceBuilderMultTranspose.hpp>
#include <sgpp/base/opencl/OCLBufferWrapperSD.hpp>
#include <sgpp/base/opencl/OCLClonedBufferMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLStretchedBufferMultiPlatform.hpp>
#include <sgpp/base/tools/QueueLoadBalancerOpenMP.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {
namespace StreamingOCLMultiPlatform {

/**
 * Kernel that provide the transposed MultiEval operation \f$v':= B v\f$ for a
 * single OpenCL
 * device.
 * This class manages the OpenCL data structures required for a OpenCL kernel
 * invocation.
 * To that end, it makes heavy use of OpenCL buffer abstraction.
 * It makes use of a queue of work packages to find out whether still device has
 * any remaining work
 * available.
 * For the creation of the device-side compute kernel code, a code generator is
 * used.
 *
 * @see base::OCLBufferWrapperSD
 * @see base::QueueLoadBalancer
 * @see SourceBuilderMultTranspose
 */
template <typename T> class KernelMultTranspose {
private:
  std::shared_ptr<base::OCLDevice> device;

  size_t dims;

  cl_int err;

  base::OCLBufferWrapperSD<T> deviceLevelTranspose;
  base::OCLBufferWrapperSD<T> deviceIndexTranspose;

  base::OCLBufferWrapperSD<T> deviceDataTranspose;
  base::OCLBufferWrapperSD<T> deviceSourceTranspose;

  base::OCLBufferWrapperSD<T> deviceResultGridTranspose;

  cl_kernel kernelMultTranspose;

  double deviceTimingMultTranspose;

  SourceBuilderMultTranspose<T> kernelSourceBuilder;
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  //    std::shared_ptr<base::OCLOperationConfiguration> parameters;
  json::Node &kernelConfiguration;

  std::shared_ptr<base::QueueLoadBalancerOpenMP> queueLoadBalancerMultTranspose;

  bool verbose;

  size_t localSize;
  size_t gridBlockSize;
  size_t scheduleSize;
  size_t totalBlockSize;

  double buildDuration;

public:
  /**
   * Constructs a new KernelMultTranspose object.
   *
   * @param device The OpenCL device this kernel instance manages
   * @param dims Dimensionality of the problem
   * @param manager The OpenCL manager to reduce OpenCL boilerplate
   * @param kernelConfiguration The configuration of this specific device
   * @param queueBalancerMultTranspose Load balance for query work from for the
   * device
   */
  KernelMultTranspose(
      std::shared_ptr<base::OCLDevice> device, size_t dims,
      std::shared_ptr<base::OCLManagerMultiPlatform> manager,
      json::Node &kernelConfiguration,
      std::shared_ptr<base::QueueLoadBalancerOpenMP> queueBalancerMultTranspose)
      : device(device), dims(dims), err(CL_SUCCESS),
        deviceLevelTranspose(device), deviceIndexTranspose(device),
        deviceDataTranspose(device), deviceSourceTranspose(device),
        deviceResultGridTranspose(device), kernelMultTranspose(nullptr),
        deviceTimingMultTranspose(0.0),
        kernelSourceBuilder(device, kernelConfiguration, dims),
        manager(manager), kernelConfiguration(kernelConfiguration),
        queueLoadBalancerMultTranspose(queueBalancerMultTranspose),
        buildDuration(0.0) {
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
    gridBlockSize =
        kernelConfiguration["KERNEL_TRANS_GRID_BLOCK_SIZE"].getUInt();
    scheduleSize = kernelConfiguration["KERNEL_SCHEDULE_SIZE"].getUInt();
    totalBlockSize = localSize * gridBlockSize;
  }

  /**
   * Destructor
   */
  ~KernelMultTranspose() {
    if (this->kernelMultTranspose != nullptr) {
      clReleaseKernel(this->kernelMultTranspose);
      this->kernelMultTranspose = nullptr;
    }
  }

  /**
   * Perform the transposed MultiEval operator \f$v':= B v\f$ with the device
   * this kernel manages.
   * Has additional, currently unused parameters to enable further MPI
   * parallelization in the
   * future.
   *
   * @param level Vector containing the d-dimensional levels of the grid, the
   * order matches the
   * index vector
   * @param index Vector containing the d-dimensional indices of the grid, the
   * order matches the
   * level vector
   * @param dataset Vector containing the d-dimensional data points
   * @param source Vector \f$v\f$
   * @param result The results vector \f$v'\f$ of the operation
   * @param start_index_data start of range of data points to work on, currently
   * not used
   * @param end_index_data end of range of data points to work on, currently not
   * used
   */
  double multTranspose(std::vector<T> &level, std::vector<T> &index,
                       std::vector<T> &dataset, std::vector<T> &source,
                       std::vector<T> &result, const size_t start_index_data,
                       const size_t end_index_data) {
    if (verbose) {
#pragma omp critical(StreamingOCLMultiPlatformKernelMultTranspose)
      {
        std::cout << "entering multTranspose, device: " << device->deviceName
                  << " (" << device->deviceId << ")" << std::endl;
      }
    }

    //    // check if there is something to do at all
    //    if (!(end_index_grid > start_index_grid && end_index_data >
    //    start_index_data)) {
    //      return 0.0;
    //    }

    if (this->kernelMultTranspose == nullptr) {
      std::chrono::time_point<std::chrono::system_clock> start, end;
      start = std::chrono::system_clock::now();

      std::string program_src = kernelSourceBuilder.generateSource();
      this->kernelMultTranspose = manager->buildKernel(
          program_src, device, kernelConfiguration, "multTransOCL");

      end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;
      this->buildDuration = elapsed_seconds.count();
    } else {
      buildDuration = 0.0;
    }

    this->deviceTimingMultTranspose = 0.0;

// for slow devices to catch up
#pragma omp barrier

    while (true) {
      size_t kernelStartData = start_index_data;
      size_t kernelEndData = end_index_data;

      // set kernel arguments
      size_t kernelStartGrid = 0;
      size_t kernelEndGrid = 0;

      bool segmentAvailable = queueLoadBalancerMultTranspose->getNextSegment(
          scheduleSize, kernelStartGrid, kernelEndGrid);
      if (!segmentAvailable) {
        break;
      }

      size_t rangeSize = kernelEndGrid - kernelStartGrid;
      size_t rangeSizeAfterBlocking =
          (kernelEndGrid / gridBlockSize) - (kernelStartGrid / gridBlockSize);

      size_t sourceSize = end_index_data - start_index_data;

      if (verbose) {
#pragma omp critical(StreamingOCLMultiPlatformKernelMultTranspose)
        {
          std::cout << "device: " << device->deviceName << " ("
                    << device->deviceId << ") "
                    << " kernel from: " << kernelStartGrid
                    << " to: " << kernelEndGrid << " -> range: " << rangeSize
                    << " (with blocking: " << rangeSizeAfterBlocking << ")"
                    << std::endl;
        }
      }

      initGridBuffersTranspose(level, index, kernelStartGrid, kernelEndGrid);
      initDatasetBuffersTranspose(dataset, source, kernelStartData,
                                  kernelEndData);
      initGridResultBuffersTranspose(kernelStartGrid, kernelEndGrid);

      clFinish(device->commandQueue);

      if (rangeSize > 0) {
        err = clSetKernelArg(kernelMultTranspose, 0, sizeof(cl_mem),
                             this->deviceLevelTranspose.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(kernelMultTranspose, 1, sizeof(cl_mem),
                             this->deviceIndexTranspose.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(kernelMultTranspose, 2, sizeof(cl_mem),
                             this->deviceDataTranspose.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(kernelMultTranspose, 3, sizeof(cl_mem),
                             this->deviceSourceTranspose.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(kernelMultTranspose, 4, sizeof(cl_mem),
                             this->deviceResultGridTranspose.getBuffer());
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err =
            clSetKernelArg(kernelMultTranspose, 5, sizeof(cl_int), &sourceSize);
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(kernelMultTranspose, 6, sizeof(cl_int),
                             &kernelStartData);
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }
        err = clSetKernelArg(kernelMultTranspose, 7, sizeof(cl_int),
                             &kernelEndData);
        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString
              << "OCL Error: Failed to create kernel arguments for device "
              << std::endl;
          throw base::operation_exception(errorString.str());
        }

        cl_event clTiming;

        err = clEnqueueNDRangeKernel(device->commandQueue, kernelMultTranspose,
                                     1, 0, &rangeSizeAfterBlocking, &localSize,
                                     0, nullptr, &clTiming);

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

        deviceResultGridTranspose.readFromBuffer();

        clFinish(device->commandQueue);
        //                std::cout << "read from buffer: " << device->deviceId
        //                << "" << std::endl;

        std::vector<T> &deviceResultGridTransposeHost =
            deviceResultGridTranspose.getHostPointer();
        for (size_t i = 0; i < rangeSize; i++) {
          result[kernelStartGrid + i] = deviceResultGridTransposeHost[i];
        }

        // determine kernel execution time
        cl_ulong startTime = 0;
        cl_ulong endTime = 0;

        err = clGetEventProfilingInfo(clTiming, CL_PROFILING_COMMAND_START,
                                      sizeof(cl_ulong), &startTime, nullptr);

        if (err != CL_SUCCESS) {
          std::stringstream errorString;
          errorString << "OCL Error: Failed to read start-time from command "
                         "queue (or crash in multTranspose)! Error code: "
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

        double time = 0.0;
        time = static_cast<double>(endTime - startTime);
        time *= 1e-9;

        if (verbose) {
#pragma omp critical(StreamingOCLMultiPlatformKernelMultTranspose)
          {
            std::cout << "device: " << device->deviceName << " ("
                      << device->deviceId << ") "
                      << "duration: " << time << std::endl;
          }
        }

        this->deviceTimingMultTranspose += time;

        clReleaseEvent(clTiming);
      }
    }
    return this->deviceTimingMultTranspose;
  }

  /**
   * @return The time it took to compile the OpenCL kernel code
   */
  double getBuildDuration() { return this->buildDuration; }

private:
  /**
   * Initializes the device buffers that belong to the grid on the device.
   * This is expensive, as it triggers host to device memory copies.
   *
   * @param level Levels vectors of the grid
   * @param index Index vectors of the grid
   * @param kernelStartGrid start of the range to be processed by the device
   * @param kernelEndGrid End of the range to be processed by the device
   */
  void initGridBuffersTranspose(std::vector<T> &level, std::vector<T> &index,
                                size_t kernelStartGrid, size_t kernelEndGrid) {
    deviceLevelTranspose.intializeTo(level, dims, kernelStartGrid,
                                     kernelEndGrid);
    deviceIndexTranspose.intializeTo(index, dims, kernelStartGrid,
                                     kernelEndGrid);
  }

  /**
   * Initializes the device buffers that belong to the dataset on the device.
   * This is expensive, as it triggers host to device memory copies.
   *
   * @param dataset The dataset to be processed
   * @param source The input vector \f$v\f$ that is transferred to the device as
   * well
   * @param kernelStartData Start of the range to be processed by the device
   * @param kernelEndData End of the range to be processed by the device
   */
  void initDatasetBuffersTranspose(std::vector<T> &dataset,
                                   std::vector<T> &source,
                                   size_t kernelStartData,
                                   size_t kernelEndData) {
    deviceDataTranspose.intializeTo(dataset, dims, kernelStartData,
                                    kernelEndData, true);
    deviceSourceTranspose.intializeTo(source, 1, kernelStartData,
                                      kernelEndData);
  }

  /**
   * Initializes the device buffers that contain the result of the multTranspose
   * operation.
   * This is expensive, as it triggers host to device memory copies.
   *
   * @param kernelStartGrid Start of the range to be processed by the device
   * @param kernelEndGrid End of the range to be processed by the device
   */
  void initGridResultBuffersTranspose(size_t kernelStartGrid,
                                      size_t kernelEndGrid) {
    size_t range = kernelEndGrid - kernelStartGrid;

    std::vector<T> zeros(range);
    for (size_t i = 0; i < range; i++) {
      zeros[i] = 0.0;
    }

    deviceResultGridTranspose.intializeTo(zeros, 1, 0, range);
  }
};

} // namespace StreamingOCLMultiPlatform
} // namespace datadriven
} // namespace sgpp

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#pragma once

#include <CL/cl.h>

#include <sgpp/base/opencl/OCLBufferWrapperSD.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

#include <limits>
#include <string>
#include <vector>

#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/KernelSourceBuilderCreateGraph.hpp>

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

/// OpenCL kernel class for operation to create a k nearest neighbor graph
template <typename T>
class KernelCreateGraph {
 private:
  /// Used opencl device
  std::shared_ptr<base::OCLDevice> device;
  size_t dims;
  size_t k;
  cl_int err;
  /// OpenCL buffer for the dataset
  base::OCLBufferWrapperSD<T> deviceData;
  /// OpenCL Buffer for the result vector
  base::OCLBufferWrapperSD<int> deviceResultData;
  cl_kernel kernel;
  /// Source builder for the opencl source code of the kernel
  sgpp::datadriven::DensityOCLMultiPlatform::SourceBuilderCreateGraph<T> kernelSourceBuilder;
  /// OpenCL manager
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  /// Stores the running time of the kernel
  double deviceTimingMult;
  bool verbose;
  /// OpenCL configuration containing the building flags
  json::Node &kernelConfiguration;
  size_t localSize;
  size_t dataBlockingSize;
  size_t scheduleSize;
  size_t totalBlockSize;
  /// Host side buffer for the dataset
  std::vector<T> &data;
  size_t unpadded_datasize;
  cl_event clTiming = nullptr;

 public:
  KernelCreateGraph(std::shared_ptr<base::OCLDevice> dev, size_t dims, size_t k,
                    std::vector<T> &data, std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                    json::Node &kernelConfiguration)
      : device(dev),
        dims(dims),
        k(k),
        err(CL_SUCCESS),
        deviceData(device),
        deviceResultData(device),
        kernel(nullptr),
        kernelSourceBuilder(kernelConfiguration, dims),
        manager(manager),
        deviceTimingMult(0.0),
        verbose(false),
        kernelConfiguration(kernelConfiguration),
        data(data) {
    this->verbose = kernelConfiguration["VERBOSE"].getBool();

    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0 &&
        kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt() < dims) {
      std::stringstream errorString;
      errorString << "OCL Error: setting \"KERNEL_DATA_STORE\" to \"register\" "
                  << "requires value of \"KERNEL_MAX_DIM_UNROLL\"";
      errorString << " to be greater than the dimension of the data set, was set to"
                  << kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt() << "(device: \""
                  << device->deviceName << "\")" << std::endl;
      throw base::operation_exception(errorString.str());
    }

    localSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    dataBlockingSize = kernelConfiguration["KERNEL_DATA_BLOCKING_SIZE"].getUInt();
    scheduleSize = kernelConfiguration["KERNEL_SCHEDULE_SIZE"].getUInt();

    if (kernelConfiguration.contains("APPROX_REG_COUNT")) {
      size_t approxRegCount = kernelConfiguration["APPROX_REG_COUNT"].getUInt();
      // Check range and whether x is a power of 2 or not
      if (approxRegCount < k || approxRegCount > localSize ||
          (approxRegCount & (approxRegCount - 1))) {
        std::stringstream errorString;
        errorString << "OCL Error: APPROX_REG_COUNT: " << approxRegCount
                    << " is not a valid size!\n"
                    << "Needs to be a power of 2, greater than k and smaller than (or equal to) the"
                    << " parameter LOCAL_SIZE: " << localSize;
        throw base::operation_exception(errorString.str());
      }
    }

    totalBlockSize = dataBlockingSize * localSize;
    unpadded_datasize = data.size();
    for (size_t i = 0; i < (localSize - data.size() % localSize) * dims; i++) data.push_back(2.0);
    deviceData.intializeTo(data, 1, 0, data.size());
  }

  ~KernelCreateGraph() {
    if (kernel != nullptr) {
      clReleaseKernel(kernel);
      this->kernel = nullptr;
    }
  }

  /// Runs the opencl kernel to find the k nearest neighbors of all datapoints in the given chunk
  void begin_graph_creation(size_t startid, size_t chunksize) {
    if (verbose) {
      std::cout << "entering graph, device: " << device->deviceName << " (" << device->deviceId
                << ")" << std::endl;
      std::cout << "k: " << k << " Dims:" << dims << std::endl;
    }
    size_t datasize = unpadded_datasize / dims;

    size_t globalworkrange[1];
    if (chunksize == 0) {
      globalworkrange[0] = unpadded_datasize / dims;
    } else {
      globalworkrange[0] = chunksize;
    }
    globalworkrange[0] = globalworkrange[0] + (localSize - globalworkrange[0] % localSize);

    // Build kernel if not already done
    if (this->kernel == nullptr) {
      if (verbose) std::cout << "generating kernel source" << std::endl;
      std::string program_src = kernelSourceBuilder.generateSource(
          dims, k, datasize,
          (unpadded_datasize / dims) + (localSize - (unpadded_datasize / dims) % localSize));
      if (verbose) std::cout << "Source: " << std::endl << program_src << std::endl;
      if (verbose) std::cout << "building kernel" << std::endl;
      this->kernel =
          manager->buildKernel(program_src, device, kernelConfiguration, "connectNeighbors");
    }

    if (!deviceResultData.isInitialized())
      deviceResultData.initializeBuffer(globalworkrange[0] * k);
    clFinish(device->commandQueue);
    this->deviceTimingMult = 0.0;

    // Set kernel arguments
    err = clSetKernelArg(this->kernel, 0, sizeof(cl_mem), this->deviceData.getBuffer());
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    err = clSetKernelArg(this->kernel, 1, sizeof(cl_mem), this->deviceResultData.getBuffer());
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    err = clSetKernelArg(this->kernel, 2, sizeof(cl_uint), &startid);
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }

    clTiming = nullptr;

    if (verbose)
      std::cout << "Starting the kernel for " << globalworkrange[0] << " items" << std::endl;
    err = clEnqueueNDRangeKernel(device->commandQueue, this->kernel, 1, 0, globalworkrange,
                                 &localSize, 0, nullptr, &clTiming);
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to enqueue kernel command! Error code: " << err
                  << std::endl;
      throw base::operation_exception(errorString.str());
    }
  }
  double finalize_graph_creation(std::vector<int> &result, size_t startid, size_t chunksize) {
    clFinish(device->commandQueue);

    if (verbose) std::cout << "Finished kernel execution" << std::endl;
    deviceResultData.readFromBuffer();
    clFinish(device->commandQueue);

    std::vector<int> &hostTemp = deviceResultData.getHostPointer();

    size_t real_count;
    if (chunksize == 0) {
      real_count = unpadded_datasize / dims;
    } else {
      real_count = chunksize;
    }
    for (size_t i = 0; i < real_count * k; i++) result[i] = hostTemp[i];
    if (verbose) std::cout << "Read data from opencl device" << std::endl;
    // determine kernel execution time
    cl_ulong startTime = 0;
    cl_ulong endTime = 0;

    err = clGetEventProfilingInfo(clTiming, CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
                                  &startTime, nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to read start-time from command queue "
                  << "(or crash in mult)! Error code: " << err << std::endl;
      throw base::operation_exception(errorString.str());
    }

    err = clGetEventProfilingInfo(clTiming, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime,
                                  nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to read end-time from command queue! "
                  << " Error code: " << err << std::endl;
      throw base::operation_exception(errorString.str());
    }

    clReleaseEvent(clTiming);

    double time = 0.0;
    time = static_cast<double>(endTime - startTime);
    time *= 1e-9;

    if (verbose) {
      std::cout << "device: " << device->deviceName << " (" << device->deviceId << ") "
                << "duration: " << time << std::endl;
    }

    this->deviceTimingMult += time;
    return 0;
  }
  /// Adds all possible building parameters to the configuration if they do not exist yet
  static void augmentDefaultParameters(sgpp::base::OCLOperationConfiguration &parameters) {
    for (std::string &platformName : parameters["PLATFORMS"].keys()) {
      json::Node &platformNode = parameters["PLATFORMS"][platformName];
      for (std::string &deviceName : platformNode["DEVICES"].keys()) {
        json::Node &deviceNode = platformNode["DEVICES"][deviceName];

        const std::string &kernelName = "connectNeighbors";

        json::Node &kernelNode = deviceNode["KERNELS"].contains(kernelName)
                                     ? deviceNode["KERNELS"][kernelName]
                                     : deviceNode["KERNELS"].addDictAttr(kernelName);

        if (kernelNode.contains("VERBOSE") == false) {
          kernelNode.addIDAttr("VERBOSE", false);
        }

        if (kernelNode.contains("LOCAL_SIZE") == false) {
          kernelNode.addIDAttr("LOCAL_SIZE", UINT64_C(128));
        }

        if (kernelNode.contains("KERNEL_USE_LOCAL_MEMORY") == false) {
          kernelNode.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
        }

        if (kernelNode.contains("KERNEL_STORE_DATA") == false) {
          kernelNode.addTextAttr("KERNEL_STORE_DATA", "array");
        }

        if (kernelNode.contains("KERNEL_MAX_DIM_UNROLL") == false) {
          kernelNode.addIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
        }

        if (kernelNode.contains("KERNEL_DATA_BLOCKING_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_DATA_BLOCKING_SIZE", UINT64_C(1));
        }

        if (kernelNode.contains("KERNEL_TRANS_GRID_BLOCKING_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_TRANS_GRID_BLOCKING_SIZE", UINT64_C(1));
        }

        if (kernelNode.contains("KERNEL_SCHEDULE_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_SCHEDULE_SIZE", UINT64_C(102400));
        }

        if (kernelNode.contains("USE_SELECT") == false) {
          kernelNode.addIDAttr("USE_SELECT", false);
        }
      }
    }
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/LinearLoadBalancerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLBufferWrapperSD.hpp>

#include <CL/cl.h>
#include <string.h>
#include <limits>
#include <string>
#include <vector>

#include "KernelSourceBuilderB.hpp"

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

template<typename T>
class KernelDensityB {
 private:
  std::shared_ptr<base::OCLDevice> device;

  size_t dims;
  size_t dataSize;
  size_t gridSize;

  cl_int err;

  base::OCLBufferWrapperSD<int> devicePoints;
  base::OCLBufferWrapperSD<T> deviceData;
  base::OCLBufferWrapperSD<T> deviceResultData;

  cl_kernel kernelB;

  sgpp::datadriven::DensityOCLMultiPlatform::SourceBuilderB<T> kernelSourceBuilder;

  std::shared_ptr<base::OCLManagerMultiPlatform> manager;

  double deviceTimingMult;

  json::Node &kernelConfiguration;


  bool verbose;

  size_t localSize;
  size_t dataBlockingSize;
  size_t scheduleSize;
  size_t totalBlockSize;

 public:
  KernelDensityB(std::shared_ptr<base::OCLDevice> dev, size_t dims,
                 std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                 json::Node &kernelConfiguration,
                 std::vector<int> &points) :
      device(dev), dims(dims), err(CL_SUCCESS), devicePoints(device),
      deviceData(device), deviceResultData(device), kernelB(nullptr),
      kernelSourceBuilder(device, kernelConfiguration, dims), manager(manager),
      deviceTimingMult(0.0), kernelConfiguration(kernelConfiguration) {
    gridSize = points.size()/(2*dims);
    this->verbose = kernelConfiguration["VERBOSE"].getBool();

    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0
        && kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt() < dims) {
      std::stringstream errorString;
      errorString
          << "OCL Error: setting \"KERNEL_DATA_STORE\" to \"register\" "
          << "requires value of \"KERNEL_MAX_DIM_UNROLL\"";
      errorString << " to be greater than the dimension of the data set, was set to"
                  << kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt() << "(device: \""
                  << device->deviceName << "\")" << std::endl;
      throw base::operation_exception(errorString.str());
    }

    localSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    dataBlockingSize = kernelConfiguration["KERNEL_DATA_BLOCKING_SIZE"].getUInt();
    scheduleSize = kernelConfiguration["KERNEL_SCHEDULE_SIZE"].getUInt();
    totalBlockSize = dataBlockingSize * localSize;

    devicePoints.intializeTo(points, 1, 0, points.size());
  }

  ~KernelDensityB() {
    if (kernelB != nullptr) {
      clReleaseKernel(kernelB);
      this->kernelB = nullptr;
    }
  }

  void resetKernel() {
  }

  double rhs(std::vector<T> &data, std::vector<T> &result,
             size_t startid = 0, size_t chunksize = -1) {
    if (verbose) {
      std::cout << "entering mult, device: " << device->deviceName << " ("
                << device->deviceId << ")" << std::endl;
    }
    this->dataSize = data.size()/dims;

    // Build kernel if not already done
    if (this->kernelB == nullptr) {
      if (verbose)
        std::cout << "generating kernel source" << std::endl;
      std::string program_src = kernelSourceBuilder.generateSource(dims, dataSize);
      if (verbose)
        std::cout << "Source: " << std::endl << program_src << std::endl;
      if (verbose)
        std::cout << "building kernel" << std::endl;
      this->kernelB = manager->buildKernel(program_src, device, kernelConfiguration, "cscheme");
    }

    // Load data into buffers if not already done
    if (!deviceData.isInitialized())
      deviceData.intializeTo(data, 1, 0, data.size());
    if (chunksize == 1) {
      std::vector<T> zeros(gridSize);
      for (size_t i = 0; i < gridSize; i++) {
        zeros[i] = 0.0;
        deviceResultData.intializeTo(zeros, 1, 0, gridSize);
      }
    } else {
      std::vector<T> zeros(chunksize);
      for (size_t i = 0; i < chunksize; i++) {
        zeros[i] = 0.0;
      }
      deviceResultData.intializeTo(zeros, 1, 0, chunksize);
    }
    clFinish(device->commandQueue);
    this->deviceTimingMult = 0.0;

    // Set kernel arguments
    err = clSetKernelArg(this->kernelB, 0, sizeof(cl_mem), this->devicePoints.getBuffer());
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    err = clSetKernelArg(this->kernelB, 1, sizeof(cl_mem), this->deviceData.getBuffer());
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    err = clSetKernelArg(this->kernelB, 2, sizeof(cl_mem), this->deviceResultData.getBuffer());
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    err = clSetKernelArg(this->kernelB, 3, sizeof(cl_uint), &startid);
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }

    cl_event clTiming = nullptr;

    // enqueue kernel
    if (verbose)
      std::cout << "Starting the kernel" << std::endl;
    size_t globalworkrange[1];
    if (chunksize == -1) {
      globalworkrange[0] = gridSize;
    } else {
      globalworkrange[0] = chunksize;
    }

    err = clEnqueueNDRangeKernel(device->commandQueue, this->kernelB, 1, 0, globalworkrange,
                                 NULL, 0, nullptr, &clTiming);
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to enqueue kernel command! Error code: " << err
                  << std::endl;
      throw base::operation_exception(errorString.str());
    }
    clFinish(device->commandQueue);

    if (verbose)
      std::cout << "Finished kernel execution" << std::endl;
    deviceResultData.readFromBuffer();
    clFinish(device->commandQueue);

    std::vector<T> &hostTemp = deviceResultData.getHostPointer();
    if (chunksize == 1) {
      for (size_t i = 0; i < gridSize; i++)
        result[i] = hostTemp[i];
    } else {
      for (size_t i = 0; i < chunksize; i++)
        result[i] = hostTemp[i];
    }

    // determine kernel execution time
    cl_ulong startTime = 0;
    cl_ulong endTime = 0;

    err = clGetEventProfilingInfo(clTiming, CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
                                  &startTime, nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString
          << "OCL Error: Failed to read start-time from command queue (or crash in mult)!"
          << " Error code: " << err << std::endl;
      throw base::operation_exception(errorString.str());
    }

    err = clGetEventProfilingInfo(clTiming, CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
                                  &endTime, nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to read end-time from command queue! Error code: "
                  << err << std::endl;
      throw base::operation_exception(errorString.str());
    }

    clReleaseEvent(clTiming);

    double time = 0.0;
    time = static_cast<double> (endTime - startTime);
    time *= 1e-9;

    if (verbose) {
      std::cout << "device: " << device->deviceName << " (" << device->deviceId << ") "
                << "duration: " << time << std::endl;
    }

    this->deviceTimingMult += time;
    return 0;
  }
  static void augmentDefaultParameters(sgpp::base::OCLOperationConfiguration &parameters){
    for (std::string &platformName : parameters["PLATFORMS"].keys()) {
      json::Node &platformNode = parameters["PLATFORMS"][platformName];
      for (std::string &deviceName : platformNode["DEVICES"].keys()) {
        json::Node &deviceNode = platformNode["DEVICES"][deviceName];

        const std::string &kernelName = "cscheme";

        json::Node &kernelNode =
            deviceNode["KERNELS"].contains(kernelName) ?
            deviceNode["KERNELS"][kernelName] :
            deviceNode["KERNELS"].addDictAttr(kernelName);

        if (kernelNode.contains("VERBOSE") == false) {
          kernelNode.addIDAttr("VERBOSE", false);
        }

        if (kernelNode.contains("LOCAL_SIZE") == false) {
          kernelNode.addIDAttr("LOCAL_SIZE", 128ul);
        }

        if (kernelNode.contains("KERNEL_USE_LOCAL_MEMORY") == false) {
          kernelNode.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
        }

        if (kernelNode.contains("KERNEL_STORE_DATA") == false) {
          kernelNode.addTextAttr("KERNEL_STORE_DATA", "array");
        }

        if (kernelNode.contains("KERNEL_MAX_DIM_UNROLL") == false) {
          kernelNode.addIDAttr("KERNEL_MAX_DIM_UNROLL", 10ul);
        }

        if (kernelNode.contains("KERNEL_DATA_BLOCKING_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_DATA_BLOCKING_SIZE", 1ul);
        }

        if (kernelNode.contains("KERNEL_TRANS_GRID_BLOCKING_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_TRANS_GRID_BLOCKING_SIZE", 1ul);
        }

        if (kernelNode.contains("KERNEL_SCHEDULE_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_SCHEDULE_SIZE", 102400ul);
        }
      }
    }
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp

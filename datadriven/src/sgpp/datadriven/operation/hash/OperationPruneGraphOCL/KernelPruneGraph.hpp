// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <CL/cl.h>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLBufferWrapperSD.hpp>

#include <limits>
#include <string>
#include <vector>

#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/KernelSourceBuilderPruneGraph.hpp>

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

/// OpenCL kernel for density based pruning of a graph
template<typename T>
class KernelPruneGraph {
 private:
  /// Used opencl device
  std::shared_ptr<base::OCLDevice> device;

  size_t dims;
  size_t gridSize;
  size_t dataSize;
  T treshold;
  size_t k;

  cl_int err;

  /// OpenCL buffer for the grid points
  base::OCLBufferWrapperSD<int> devicePoints;
  /// OpenCL buffer for the alpha vector used to calculate the density
  base::OCLBufferWrapperSD<T> deviceAlpha;
  /// OpenCL buffer for the dataset
  base::OCLBufferWrapperSD<T> deviceData;
  /// OpenCL buffer for the graph
  base::OCLBufferWrapperSD<int> deviceGraph;

  cl_kernel kernel;
  /// Source builder for the opencl source code of the kernel
  sgpp::datadriven::DensityOCLMultiPlatform::SourceBuilderPruneGraph<T> kernelSourceBuilder;
  /// OpenCL manager
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  /// Stores the running time of the kernel
  double deviceTimingMult;
  /// OpenCL configuration containing the building flags
  json::Node &kernelConfiguration;


  bool verbose;

  size_t localSize;
  size_t dataBlockingSize;
  size_t scheduleSize;
  size_t totalBlockSize;

 public:
  KernelPruneGraph(std::shared_ptr<base::OCLDevice> dev, size_t dims, T treshold, size_t k,
                   std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                   json::Node &kernelConfiguration,
                   std::vector<int> &pointsVector, std::vector<T> &alphaVector,
                   std::vector<T> &dataVector) :
      device(dev), dims(dims), treshold(treshold), k(k), err(CL_SUCCESS), devicePoints(device),
      deviceAlpha(device),
      deviceData(device), deviceGraph(device), kernel(nullptr),
      kernelSourceBuilder(kernelConfiguration, dims), manager(manager),
      deviceTimingMult(0.0),
      kernelConfiguration(kernelConfiguration) {
    this->verbose = kernelConfiguration["VERBOSE"].getBool();
    gridSize = pointsVector.size()/(2*dims);

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

    devicePoints.intializeTo(pointsVector, 1, 0, gridSize*dims*2);
    deviceAlpha.intializeTo(alphaVector, 1, 0, gridSize);
    deviceData.intializeTo(dataVector, 1, 0, dataVector.size());
    clFinish(device->commandQueue);
  }

  ~KernelPruneGraph() {
    if (kernel != nullptr) {
      clReleaseKernel(kernel);
      this->kernel = nullptr;
    }
  }

  /// Deletes nodes and edges in areas of low density
  double prune_graph(std::vector<int> &graph, size_t startid, size_t chunksize) {
    if (verbose)
      std::cout << "entering prune graph, device: " << device->deviceName
                << " (" << device->deviceId << ")"
                << std::endl;

    // Build kernel if not already done
    if (this->kernel == nullptr) {
      if (verbose)
        std::cout << "generating kernel source" << std::endl;
      std::string program_src = kernelSourceBuilder.generateSource(dims, gridSize,
                                                                   k, treshold);
      if (verbose)
        std::cout << "Source: " << std::endl << program_src << std::endl;
      if (verbose)
        std::cout << "building kernel" << std::endl;
      this->kernel = manager->buildKernel(program_src, device, kernelConfiguration, "removeEdges");
    }

    deviceGraph.intializeTo(graph, 1, 0, graph.size());
    clFinish(device->commandQueue);
    this->deviceTimingMult = 0.0;

    // Set kernel arguments
    err = clSetKernelArg(this->kernel, 0, sizeof(cl_mem), this->deviceGraph.getBuffer());
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    err = clSetKernelArg(this->kernel, 1, sizeof(cl_mem), this->devicePoints.getBuffer());
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    err = clSetKernelArg(this->kernel, 2, sizeof(cl_mem), this->deviceData.getBuffer());
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    err = clSetKernelArg(this->kernel, 3, sizeof(cl_mem), this->deviceAlpha.getBuffer());
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    err = clSetKernelArg(this->kernel, 4, sizeof(cl_uint), &startid);
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
    if (chunksize == 0) {
      globalworkrange[0] = graph.size()/k;
    } else {
      globalworkrange[0] = chunksize;
    }
    err = clEnqueueNDRangeKernel(device->commandQueue, this->kernel, 1, 0, globalworkrange,
                                 nullptr, 0, nullptr, &clTiming);
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to enqueue kernel command! Error code: " << err
                  << std::endl;
      throw base::operation_exception(errorString.str());
    }
    clFinish(device->commandQueue);

    if (verbose)
      std::cout << "Finished kernel execution" << std::endl;
    deviceGraph.readFromBuffer();
    clFinish(device->commandQueue);

    std::vector<int> &hostTemp = deviceGraph.getHostPointer();
    for (size_t i = 0; i < graph.size(); i++)
      graph[i] = hostTemp[i];
    // determine kernel execution time
    cl_ulong startTime = 0;
    cl_ulong endTime = 0;

    err = clGetEventProfilingInfo(clTiming, CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
                                  &startTime, nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString
          << "OCL Error: Failed to read start-time from command queue "
          << "(or crash in prune)! Error code: " << err << std::endl;
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

    if (verbose)
      std::cout << "device: " << device->deviceName << " (" << device->deviceId << ") "
                << "duration: " << time << std::endl;

    this->deviceTimingMult += time;
    return 0;
  }
  /// Adds all possible building parameters to the configuration if they do not exist yet
  static void augmentDefaultParameters(sgpp::base::OCLOperationConfiguration &parameters) {
    for (std::string &platformName : parameters["PLATFORMS"].keys()) {
      json::Node &platformNode = parameters["PLATFORMS"][platformName];
      for (std::string &deviceName : platformNode["DEVICES"].keys()) {
        json::Node &deviceNode = platformNode["DEVICES"][deviceName];

        const std::string &kernelName = "removeEdges";

        json::Node &kernelNode =
            deviceNode["KERNELS"].contains(kernelName) ?
            deviceNode["KERNELS"][kernelName] :
            deviceNode["KERNELS"].addDictAttr(kernelName);

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
      }
    }
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp

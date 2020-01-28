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

#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/KernelSourceBuilderMult.hpp>

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

/// Class for the OpenCL density matrix vector multiplication
template <typename T>
class KernelDensityMult {
 private:
  /// Used opencl device
  std::shared_ptr<base::OCLDevice> device;

  size_t dims;
  size_t gridSize;
  T lambda;

  cl_int err;

  /// Buffer for the grid points
  base::OCLBufferWrapperSD<int> devicePoints;
  /// Buffer for the alpha vector
  base::OCLBufferWrapperSD<T> deviceAlpha;
  /// Buffer for the result vector
  base::OCLBufferWrapperSD<T> deviceResultData;
  /// Buffer for the preprocessed levels (if used)
  base::OCLBufferWrapperSD<T> deviceLevels;
  /// Buffer for the preprocessed inversed levels (if used)
  base::OCLBufferWrapperSD<T> deviceDivisors;
  /// Buffer for the preprocessed h^{-1} (if used)
  base::OCLBufferWrapperSD<int> devicehInverse;
  /// Buffer for the preprocessed h (if used)
  base::OCLBufferWrapperSD<T> devicehs;
  /// Buffer for the preprocessed h (if used)
  base::OCLBufferWrapperSD<T> devicePositions;

  cl_kernel kernelMult;
  /// Source builder for the opencl source code of the kernel
  SourceBuilderMult<T> kernelSourceBuilder;
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

  /// Use a cache for the 2^l values? Configuration parameter is USE_LEVEL_CACHE
  bool use_level_cache;
  /// Use a calculation scheme with less operations but more branching?
  bool use_less;
  /// Use ternary operator for branching? Configuration parameter is USE_LESS_OPERATIONS
  bool do_not_use_ternary;
  /// Use preprocessed grid positions? Configuration parameter is PREPROCESSED_POSITIONS
  bool preprocess_positions;

  cl_event clTiming;

 public:
  KernelDensityMult(std::shared_ptr<base::OCLDevice> dev, size_t dims,
                    std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                    json::Node &kernelConfiguration, std::vector<int> &points, T lambda)
      : device(dev),
        dims(dims),
        lambda(lambda),
        err(CL_SUCCESS),
        devicePoints(device),
        deviceAlpha(device),
        deviceResultData(device),
        deviceLevels(device),
        deviceDivisors(device),
        devicehInverse(device),
        devicehs(device),
        devicePositions(device),
        kernelMult(nullptr),
        kernelSourceBuilder(kernelConfiguration),
        manager(manager),
        deviceTimingMult(0.0),
        kernelConfiguration(kernelConfiguration),
        use_level_cache(false),
        use_less(false),
        do_not_use_ternary(false),
        preprocess_positions(false) {
    this->verbose = use_level_cache = kernelConfiguration["VERBOSE"].getBool();
    gridSize = points.size() / (2 * dims);
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

    // Get local size and push grid points into opencl buffer
    localSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    dataBlockingSize = kernelConfiguration["KERNEL_DATA_BLOCKING_SIZE"].getUInt();
    scheduleSize = kernelConfiguration["KERNEL_SCHEDULE_SIZE"].getUInt();
    totalBlockSize = dataBlockingSize * localSize;
    for (size_t i = 0; i < (localSize - (gridSize % localSize)) * 2 * dims; i++) {
      points.push_back(0);
    }
    devicePoints.intializeTo(points, 1, 0, points.size());

    // Check whether to calculate h_n on the host side (true) or in the opencl kernel (false)
    if (kernelConfiguration.contains("USE_LEVEL_CACHE")) {
      use_level_cache = kernelConfiguration["USE_LEVEL_CACHE"].getBool();
      // Do not use level cache if preprocessed positions are enabled
      if (kernelConfiguration.contains("PREPROCESS_POSITIONS")) {
        if (kernelConfiguration["PREPROCESS_POSITIONS"].getBool()) use_level_cache = false;
      }
    }
    if (use_level_cache) {
      // Find biggest level
      int maxlevel = 1;
      for (size_t gridpoint = 0; gridpoint < gridSize; gridpoint++) {
        for (size_t d = 0; d < dims; ++d) {
          if (points[gridpoint * dims + d * 2 + 1] > maxlevel) {
            maxlevel = points[gridpoint * dims + d * 2 + 1];
          }
        }
      }
      // Store h
      std::vector<T> hs;
      hs.push_back(0.0);
      for (int l = 1; l <= maxlevel; ++l) {
        hs.push_back(static_cast<T>(1.0 / (1 << l)));
      }
      deviceLevels.intializeTo(hs, 1, 0, maxlevel + 1);
    }

    // Check whether to use more operations to avoid branching on integrating base functions
    // with same level
    if (kernelConfiguration.contains("USE_LESS_OPERATIONS"))
      use_less = kernelConfiguration["USE_LESS_OPERATIONS"].getBool();
    // Checking whether to use ternary operator for integrating base functions with
    // the same level
    if (kernelConfiguration.contains("DO_NOT_USE_TERNARY"))
      do_not_use_ternary = kernelConfiguration["DO_NOT_USE_TERNARY"].getBool();
    if (do_not_use_ternary) {
      std::vector<T> divisors;
      T divisor = 1.0;
      for (size_t i = 0; i < dims + 1; ++i) {
        divisors.push_back(divisor);
        divisor /= static_cast<T>(3.0);
      }
      deviceDivisors.intializeTo(divisors, 1, 0, dims + 1);
    }

    // Check whether to use prepocessed positions or not
    if (kernelConfiguration.contains("PREPROCESS_POSITIONS")) {
      preprocess_positions = kernelConfiguration["PREPROCESS_POSITIONS"].getBool();
    }
    if (preprocess_positions) {
      std::vector<T> positions;
      std::vector<T> hs;
      std::vector<int> hs_inverse;
      for (size_t i = 0; i < points.size() / (2 * dims); ++i) {
        for (size_t dim = 0; dim < dims; ++dim) {
          hs_inverse.push_back(1 << points[i * 2 * dims + 2 * dim + 1]);
          hs.push_back(static_cast<T>(1.0 / (1 << points[i * 2 * dims + 2 * dim + 1])));
          positions.push_back(static_cast<T>(points[i * 2 * dims + 2 * dim]) /
                              static_cast<T>(1 << points[i * 2 * dims + 2 * dim + 1]));
        }
      }
      devicehInverse.intializeTo(hs_inverse, 1, 0, points.size() / 2);
      devicehs.intializeTo(hs, 1, 0, points.size() / 2);
      devicePositions.intializeTo(positions, 1, 0, points.size() / 2);
    }

    // Finish writing all buffers
    clFinish(device->commandQueue);
  }

  ~KernelDensityMult() {
    // Release kernel
    if (kernelMult != nullptr) {
      clReleaseKernel(kernelMult);
      this->kernelMult = nullptr;
    }
  }

  void initialize_alpha_buffer(std::vector<T> &alpha) {
    // Load data into buffers if not already done
    for (size_t i = 0; i < localSize - gridSize % localSize; i++) alpha.push_back(0.0);
    deviceAlpha.intializeTo(alpha, 1, 0, alpha.size());
  }

  /// Executes part of one matrix-vector density multiplication from startid to startid + chunksize
  void start_mult(size_t startid, size_t chunksize) {
    if (verbose) {
      std::cout << "entering mult, device: " << device->deviceName << " (" << device->deviceId
                << ")" << std::endl;
    }

    // Calculate workrange for current multiplication
    size_t globalworkrange[1];
    if (chunksize == 0) {
      globalworkrange[0] = gridSize / dataBlockingSize;
    } else {
      globalworkrange[0] = chunksize / dataBlockingSize;
    }
    if (globalworkrange[0] % localSize != 0)
      globalworkrange[0] = globalworkrange[0] + (localSize - globalworkrange[0] % localSize);

    // Generate OpenCL source and build kernel if not already done
    if (this->kernelMult == nullptr) {
      if (verbose) std::cout << "generating kernel source" << std::endl;
      std::string program_src =
          kernelSourceBuilder.generateSource(dims, gridSize + (localSize - gridSize % localSize));
      if (verbose) std::cout << "Source: " << std::endl << program_src << std::endl;
      if (verbose) std::cout << "building kernel" << std::endl;
      this->kernelMult =
          manager->buildKernel(program_src, device, kernelConfiguration, "multdensity");
    }

    if (!deviceResultData.isInitialized()) {
      if (chunksize == 0)
        deviceResultData.initializeBuffer(gridSize + localSize - gridSize % localSize);
      else
        deviceResultData.initializeBuffer(globalworkrange[0] * dataBlockingSize);
      this->deviceTimingMult = 0.0;
      clFinish(device->commandQueue);
    }
    // Set mandatory kernel arguments
    int argument_counter = 0;
    if (!preprocess_positions) {
      err = clSetKernelArg(this->kernelMult, argument_counter, sizeof(cl_mem),
                           this->devicePoints.getBuffer());
      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to create kernel arguments (argument " << argument_counter
                    << ") for device " << std::endl;
        throw base::operation_exception(errorString.str());
      }
      argument_counter++;
    } else {
      err = clSetKernelArg(this->kernelMult, argument_counter, sizeof(cl_mem),
                           devicehInverse.getBuffer());
      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to create kernel arguments (argument " << argument_counter
                    << ") for device " << std::endl;
        throw base::operation_exception(errorString.str());
      }
      argument_counter++;
      err =
          clSetKernelArg(this->kernelMult, argument_counter, sizeof(cl_mem), devicehs.getBuffer());
      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to create kernel arguments (argument " << argument_counter
                    << ") for device " << std::endl;
        throw base::operation_exception(errorString.str());
      }
      argument_counter++;
      err = clSetKernelArg(this->kernelMult, argument_counter, sizeof(cl_mem),
                           devicePositions.getBuffer());
      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to create kernel arguments (argument " << argument_counter
                    << ") for device " << std::endl;
        throw base::operation_exception(errorString.str());
      }
      argument_counter++;
    }
    err = clSetKernelArg(this->kernelMult, argument_counter, sizeof(cl_mem),
                         this->deviceAlpha.getBuffer());
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments (argument " << argument_counter
                  << ") for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    argument_counter++;
    err = clSetKernelArg(this->kernelMult, argument_counter, sizeof(cl_mem),
                         this->deviceResultData.getBuffer());
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments (argument " << argument_counter
                  << ") for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    argument_counter++;
    if (std::is_same<T, float>::value) {
      err = clSetKernelArg(this->kernelMult, argument_counter, sizeof(cl_float), &lambda);
    } else {
      err = clSetKernelArg(this->kernelMult, argument_counter, sizeof(cl_double), &lambda);
    }
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments (argument " << argument_counter
                  << ") for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    argument_counter++;
    err = clSetKernelArg(this->kernelMult, argument_counter, sizeof(cl_uint), &startid);
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel arguments (argument " << argument_counter
                  << ") for device " << std::endl;
      throw base::operation_exception(errorString.str());
    }
    argument_counter++;

    // Set optional kernel arguments
    if (use_level_cache) {
      err = clSetKernelArg(this->kernelMult, argument_counter, sizeof(cl_mem),
                           this->deviceLevels.getBuffer());
      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to create kernel arguments (argument " << argument_counter
                    << ") for device " << std::endl;
        throw base::operation_exception(errorString.str());
      }
      argument_counter++;
    }
    if (do_not_use_ternary) {
      err = clSetKernelArg(this->kernelMult, argument_counter, sizeof(cl_mem),
                           this->deviceDivisors.getBuffer());
      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to create kernel arguments (argument " << argument_counter
                    << ") for device " << std::endl;
        throw base::operation_exception(errorString.str());
      }
    }
    clTiming = nullptr;
    // enqueue kernel
    if (verbose)
      std::cout << "Starting the kernel with " << globalworkrange[0] << "items" << std::endl;
    err = clEnqueueNDRangeKernel(device->commandQueue, this->kernelMult, 1, 0, globalworkrange,
                                 &localSize, 0, nullptr, &clTiming);
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to enqueue kernel command! Error code: " << err
                  << std::endl;
      throw base::operation_exception(errorString.str());
    }
    return;
  }

  double finish_mult(std::vector<T> &result, int startid, int chunksize) {
    clFinish(device->commandQueue);

    if (verbose) std::cerr << "Finished kernel execution" << std::endl;
    deviceResultData.readFromBuffer();
    clFinish(device->commandQueue);

    std::vector<T> &hostTemp = deviceResultData.getHostPointer();
    if (chunksize == 0) {
      for (size_t i = 0; i < gridSize; i++) result[i] = hostTemp[i];
    } else {
      for (int i = 0; i < chunksize; i++) result[i] = hostTemp[i];
    }
    // determine kernel execution time
    cl_ulong startTime = 0;
    cl_ulong endTime = 0;

    err = clGetEventProfilingInfo(clTiming, CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
                                  &startTime, nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to read start-time from command queue (or crash in mult)!"
                  << " Error code: " << err << std::endl;
      throw base::operation_exception(errorString.str());
    }

    err = clGetEventProfilingInfo(clTiming, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime,
                                  nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to read end-time from command queue! Error code: " << err
                  << std::endl;
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

        const std::string &kernelName = "multdensity";

        json::Node &kernelNode = deviceNode["KERNELS"].contains(kernelName)
                                     ? deviceNode["KERNELS"][kernelName]
                                     : deviceNode["KERNELS"].addDictAttr(kernelName);

        if (kernelNode.contains("REUSE_SOURCE") == false) {
          kernelNode.addIDAttr("REUSE_SOURCE", false);
        }

        if (kernelNode.contains("WRITE_SOURCE") == false) {
          kernelNode.addIDAttr("WRITE_SOURCE", false);
        }

        if (kernelNode.contains("VERBOSE") == false) {
          kernelNode.addIDAttr("VERBOSE", false);
        }

        if (kernelNode.contains("LOCAL_SIZE") == false) {
          kernelNode.addIDAttr("LOCAL_SIZE", UINT64_C(128));
        }

        if (kernelNode.contains("KERNEL_USE_LOCAL_MEMORY") == false) {
          kernelNode.addIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
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
        if (kernelNode.contains("USE_LEVEL_CACHE") == false) {
          kernelNode.addIDAttr("USE_LEVEL_CACHE", false);
        }
        if (kernelNode.contains("USE_LESS_OPERATIONS") == false) {
          kernelNode.addIDAttr("USE_LESS_OPERATIONS", true);
        }

        if (kernelNode.contains("DO_NOT_USE_TERNARY") == false) {
          kernelNode.addIDAttr("DO_NOT_USE_TERNARY", false);
        }

        if (kernelNode.contains("USE_IMPLICIT") == false) {
          kernelNode.addIDAttr("USE_IMPLICIT", false);
        }

        if (kernelNode.contains("USE_FABS") == false) {
          kernelNode.addIDAttr("USE_FABS", false);
        }
        if (kernelNode.contains("PREPROCESS_POSITIONS") == false) {
          kernelNode.addIDAttr("PREPROCESS_POSITIONS", false);
        }
      }
    }
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp

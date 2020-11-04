// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

// define required for clCreateCommandQueue on platforms that don't support
// OCL2.0 yet
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#include <CL/cl.h>

#include <map>
#include <vector>
#include <string>

#include <sgpp/base/opencl/OCLPlatformWrapper.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLDevice.hpp>

namespace sgpp {
namespace base {

class OCLManagerMultiPlatform {
 public:
  std::shared_ptr<base::OCLOperationConfiguration> parameters;

  std::vector<OCLPlatformWrapper> platforms;

  // linear device list
  std::vector<std::shared_ptr<OCLDevice>> devices;

  cl_uint overallDeviceCount;  // devices over all platforms

  bool verbose;

 public:
  explicit OCLManagerMultiPlatform(bool verbose = false);

  explicit OCLManagerMultiPlatform(std::shared_ptr<base::OCLOperationConfiguration> parameters);

  ~OCLManagerMultiPlatform();

  /**
   * @brief buildKernel builds the program that is represented by @a program_src
   * and creates @a num_devices kernel objects
   * that are stored into the map @a kernels
   *
   * @param program_src the source of the program to compile
   * @param kernel_name name of the kernel function (in program_src) to create
   * the kernel for
   * @param kernels the resulting kernels are put into this map, one for each platform
   */
  void buildKernel(const std::string &program_src, const std::string &kernel_name,
                   std::map<cl_platform_id, std::vector<cl_kernel>> &kernels);

  cl_kernel buildKernel(const std::string &source, std::shared_ptr<OCLDevice> device,
                        json::Node &kernelConfiguration, const std::string &kernelName);

  void configure(base::OCLOperationConfiguration &configuration, bool useConfiguration = false);

  void configurePlatform(cl_platform_id platformId, base::OCLOperationConfiguration &configuration,
                         bool useConfiguration);

  void configureDevice(cl_device_id deviceId, json::Node &devicesNode,
                       std::vector<cl_device_id> &filteredDeviceIds,
                       std::vector<std::string> &filteredDeviceNames,
                       std::map<std::string, size_t> &countLimitMap, bool useConfiguration);

  std::shared_ptr<base::OCLOperationConfiguration> getConfiguration();

  std::vector<std::shared_ptr<OCLDevice>> &getDevices();
};
}  // namespace base
}  // namespace sgpp

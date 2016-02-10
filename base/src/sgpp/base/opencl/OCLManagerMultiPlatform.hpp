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

#include <sgpp/base/opencl/OCLPlatformWrapper.hpp>
#include "OCLOperationConfiguration.hpp"
#include "OCLDevice.hpp"

namespace SGPP {
namespace base {

class OCLManagerMultiPlatform {
 public:
  std::shared_ptr<base::OCLOperationConfiguration> parameters;

  std::vector<OCLPlatformWrapper> platforms;

  // linear device list
  std::vector<std::shared_ptr<OCLDevice>> devices;
  //        cl_uint platformCount;
  //
  //        std::map<cl_platform_id, size_t> platformDeviceCount; //devices on
  //        the individual platforms
  //        cl_platform_id *platformIds;
  //
  //        std::map<cl_platform_id, cl_device_id *> platformDeviceIds; //
  //        device ids over all platforms
  cl_uint overallDeviceCount;  // devices over all platforms

  //        //platforms -> (deviceId -> command_queue)
  //        std::map<cl_uint, std::map<cl_uint, cl_command_queue> >
  //        commandQueues;
  //        std::vector<cl_context> platformContext;
  bool verbose;

 public:
  OCLManagerMultiPlatform(bool verbose = false);

  OCLManagerMultiPlatform(
      std::shared_ptr<base::OCLOperationConfiguration> parameters);

  ~OCLManagerMultiPlatform();

  /**
   * @brief buildKernel builds the program that is represented by @a program_src
   * and creates @a num_devices kernel objects
   * that are stored into the array @a kernel (must be already allocated with at
   * least @a num_devices )
   *
   * @param program_src the source of the program to compile
   * @param kernel_name name of the kernel function (in program_src) to create
   * the kernel for
   * @param context OpenCL context
   * @param num_devices number of OpenCL devices
   * @param device_ids array with device ids, necessary for displaying build
   * info
   * @param kernel already allocated array: the resulting kernels are put into
   * this array, one for each device (=> at least num_devices entries)
   * @return
   */
  void buildKernel(const std::string &program_src,
                   const std::string &kernel_name,
                   std::map<cl_platform_id, std::vector<cl_kernel>> &kernels);

  cl_kernel buildKernel(const std::string &source,
                        std::shared_ptr<OCLDevice> device,
                        const std::string &kernelName);

  //    void setPlatformIDs();
  //
  //    void printPlatformsInfo();
  //
  //    void setupPlatforms();
  //
  //    void setTotalDeviceCount();
  //
  //    void setDeviceType();
  //
  //    void setupDeviceIDs();

  void configure(base::OCLOperationConfiguration &configuration,
                 bool useConfiguration = false);

  void configurePlatform(cl_platform_id platformId,
                         base::OCLOperationConfiguration &configuration,
                         bool useConfiguration);

  void configureDevice(cl_device_id deviceId, json::Node &devicesNode,
                       std::vector<cl_device_id> &filteredDeviceIds,
                       std::vector<std::string> &filteredDeviceNames,
                       std::map<std::string, size_t> &countLimitMap,
                       bool useConfiguration);

  std::shared_ptr<base::OCLOperationConfiguration> getConfiguration();

  std::vector<std::shared_ptr<OCLDevice>> &getDevices();
};
}
}

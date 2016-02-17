// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <vector>

#include "sgpp/base/opencl/OCLManagerMultiPlatform.hpp"
#include "sgpp/base/exception/operation_exception.hpp"

namespace SGPP {
namespace base {

OCLManagerMultiPlatform::OCLManagerMultiPlatform(bool verbose) {
  parameters = std::make_shared<OCLOperationConfiguration>();
  parameters->replaceIDAttr("VERBOSE", verbose);
  parameters->replaceIDAttr("OCL_MANAGER_VERBOSE", false);
  parameters->replaceIDAttr("SHOW_BUILD_LOG", false);
  parameters->replaceDictAttr("PLATFORMS");
  parameters->replaceIDAttr("LOAD_BALANCING_VERBOSE", false);

  // TODO(pfandedd): is it a good idea to set default precision here?
  parameters->replaceTextAttr("INTERNAL_PRECISION", "double");

  this->verbose = verbose;
  overallDeviceCount = 0;

  this->configure(*parameters, false);
}

OCLManagerMultiPlatform::OCLManagerMultiPlatform(
    std::shared_ptr<base::OCLOperationConfiguration> parameters)
    : parameters(parameters) {
  if (!parameters->contains("VERBOSE")) {
    parameters->addIDAttr("VERBOSE", false);
  }
  if (!parameters->contains("OCL_MANAGER_VERBOSE")) {
    parameters->replaceIDAttr("OCL_MANAGER_VERBOSE", false);
  }
  if (!parameters->contains("SHOW_BUILD_LOG")) {
    parameters->replaceIDAttr("SHOW_BUILD_LOG", false);
  }
  if (!parameters->contains("PLATFORMS")) {
    parameters->replaceListAttr("PLATFORMS");
  }
  if (!parameters->contains("LOAD_BALANCING_VERBOSE")) {
    parameters->replaceIDAttr("LOAD_BALANCING_VERBOSE", false);
  }

  // TODO(pfandedd): is it a good idea to set default precision here?
  if (!parameters->contains("INTERNAL_PRECISION")) {
    parameters->replaceTextAttr("INTERNAL_PRECISION", "double");
  }

  this->verbose = (*parameters)["VERBOSE"].getBool();
  this->overallDeviceCount = 0;

  this->configure(*parameters, true);
}

OCLManagerMultiPlatform::~OCLManagerMultiPlatform() {
  cl_int err;
  for (OCLPlatformWrapper platform : this->platforms) {
    for (size_t i = 0; i < platform.deviceIds.size(); i++) {
      err = clReleaseCommandQueue(platform.commandQueues[i]);
      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Could not release command queue! Error Code: " << err
                    << std::endl;
        throw SGPP::base::operation_exception(errorString.str());
      }
    }
    cl_int err = clReleaseContext(platform.context);
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Could not release context! Error Code: " << err << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }
  }
}

void OCLManagerMultiPlatform::buildKernel(
    const std::string &program_src, const std::string &kernel_name,
    std::map<cl_platform_id, std::vector<cl_kernel>> &kernels) {
  cl_int err;

  for (OCLPlatformWrapper &platform : this->platforms) {
    // setting the program
    const char *kernel_src = program_src.c_str();
    cl_program program = clCreateProgramWithSource(platform.context, 1, &kernel_src, NULL, &err);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create program! Error Code: " << err << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }

    std::string build_opts;

    if (!(*parameters).contains("ENABLE_OPTIMIZATIONS") ||
        (*parameters)["ENABLE_OPTIMIZATIONS"].getBool()) {
      // TODO(pfandedd): user should be able to change
      std::string optimizationFlags = "";
      if ((*parameters).contains("OPTIMIZATION_FLAGS")) {
        optimizationFlags = (*parameters)["OPTIMIZATION_FLAGS"].get();
      }
      build_opts = optimizationFlags;  // -O5  -cl-mad-enable -cl-denorms-are-zero
                                       // -cl-no-signed-zeros
                                       // -cl-unsafe-math-optimizations
                                       // -cl-finite-math-only -cl-fast-relaxed-math
    } else {
      build_opts = "-cl-opt-disable";  // -g
    }

    // compiling the program
    err = clBuildProgram(program, 0, NULL, build_opts.c_str(), NULL, NULL);

    if (err != CL_SUCCESS) {
      // get the build log
      size_t len;
      clGetProgramBuildInfo(program, platform.deviceIds[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
      std::string buffer(len, '\0');
      clGetProgramBuildInfo(program, platform.deviceIds[0], CL_PROGRAM_BUILD_LOG, len, &buffer[0],
                            NULL);
      buffer = buffer.substr(0, buffer.find('\0'));

      if (verbose) {
        std::cout << "--- Build Log ---" << std::endl
                  << buffer << std::endl;
      }

      std::stringstream errorString;
      errorString << "OCL Error: OpenCL build error. Error code: " << err << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }

    for (size_t i = 0; i < platform.deviceIds.size(); i++) {
      // creating the kernel
      cl_kernel kernel = clCreateKernel(program, kernel_name.c_str(), &err);
      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to create kernel! Error code: " << err << std::endl;
        throw SGPP::base::operation_exception(errorString.str());
      }

      kernels[platform.platformId].push_back(kernel);
    }

    if (program) {
      clReleaseProgram(program);
    }
  }
}

cl_kernel OCLManagerMultiPlatform::buildKernel(const std::string &source,
                                               std::shared_ptr<OCLDevice> device,
                                               json::Node &kernelConfiguration,
                                               const std::string &kernelName) {
  cl_int err;

  // setting the program
  const char *kernel_src = source.c_str();
  cl_program program = clCreateProgramWithSource(device->context, 1, &kernel_src, NULL, &err);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Failed to create program! Error Code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  std::string build_opts;
  if (!kernelConfiguration.contains("ENABLE_OPTIMIZATIONS") ||
      kernelConfiguration["ENABLE_OPTIMIZATIONS"].getBool()) {
    // TODO(pfandedd): user should be able to change
    std::string optimizationFlags = "";
    if (kernelConfiguration.contains("OPTIMIZATION_FLAGS")) {
      optimizationFlags = kernelConfiguration["OPTIMIZATION_FLAGS"].get();
    }
    build_opts = optimizationFlags;  // -O5  -cl-mad-enable -cl-denorms-are-zero
    // -cl-no-signed-zeros -cl-unsafe-math-optimizations
    // -cl-finite-math-only -cl-fast-relaxed-math
  } else {
    build_opts = "-cl-opt-disable";  // -g
  }

  // compiling the program
  err = clBuildProgram(program, 0, NULL, build_opts.c_str(), NULL, NULL);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: OpenCL build error. Error code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  // get the build log
  size_t len;
  clGetProgramBuildInfo(program, device->deviceId, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
  std::string buffer(len, '\0');
  clGetProgramBuildInfo(program, device->deviceId, CL_PROGRAM_BUILD_LOG, len, &buffer[0], NULL);
  buffer = buffer.substr(0, buffer.find('\0'));

  if (verbose) {
    std::cout << "--- Build Log ---" << std::endl
              << buffer << std::endl;
  }

  cl_kernel kernel = clCreateKernel(program, kernelName.c_str(), &err);
  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Failed to create kernel! Error code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  if (program) {
    clReleaseProgram(program);
  }
  return kernel;
}

std::shared_ptr<base::OCLOperationConfiguration> OCLManagerMultiPlatform::getConfiguration() {
  return this->parameters;
}

void OCLManagerMultiPlatform::configure(base::OCLOperationConfiguration &configuration,
                                        bool useConfiguration) {
  cl_int err;

  // determine number of available OpenCL platforms
  cl_uint platformCount;
  err = clGetPlatformIDs(0, nullptr, &platformCount);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get number of OpenCL platforms. Error Code: " << err
                << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  if (verbose) {
    std::cout << "OCL Info: " << platformCount << " OpenCL platforms have been found" << std::endl;
  }

  // get available platforms
  std::vector<cl_platform_id> platformIds(platformCount);
  err = clGetPlatformIDs(platformCount, platformIds.data(), nullptr);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get Platform ID. Error Code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  for (size_t i = 0; i < platformCount; i++) {
    this->configurePlatform(platformIds[i], *parameters, useConfiguration);
  }
}

void OCLManagerMultiPlatform::configurePlatform(cl_platform_id platformId,
                                                base::OCLOperationConfiguration &configuration,
                                                bool useConfiguration) {
  cl_int err;

  char platformName[128] = {0};
  err = clGetPlatformInfo(platformId, CL_PLATFORM_NAME, 128 * sizeof(char), platformName, nullptr);

  if (CL_SUCCESS != err) {
    std::stringstream errorString;
    errorString << "OCL Error: can't get platform name!" << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  } else {
    if (verbose) {
      std::cout << "OCL Info: detected platform, name: \"" << platformName << "\"" << std::endl;
    }
  }

  if (verbose) {
    char vendor_name[128] = {0};
    err =
        clGetPlatformInfo(platformId, CL_PLATFORM_VENDOR, 128 * sizeof(char), vendor_name, nullptr);

    if (CL_SUCCESS != err) {
      std::stringstream errorString;
      errorString << "OCL Error: Can't get platform vendor!" << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    } else {
      std::cout << "OCL Info: detected platform vendor name: \"" << vendor_name << "\""
                << std::endl;
    }
  }

  if (useConfiguration) {
    if (!(*parameters)["PLATFORMS"].contains(platformName)) {
      return;
    }
  } else {
    // creating new configuration
    json::Node &platformNode = (*parameters)["PLATFORMS"].addDictAttr(platformName);
    platformNode.addDictAttr("DEVICES");
  }

  if (verbose) {
    std::cout << "OCL Info: using platform, name: \"" << platformName << "\"" << std::endl;
  }

  json::Node &devicesNode = (*parameters)["PLATFORMS"][platformName]["DEVICES"];

  uint32_t currentPlatformDevices;
  // get the number of devices
  err = clGetDeviceIDs(platformId, CL_DEVICE_TYPE_ALL, 0, nullptr, &currentPlatformDevices);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get device count. Error Code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  std::vector<cl_device_id> deviceIds(currentPlatformDevices);
  err = clGetDeviceIDs(platformId, CL_DEVICE_TYPE_ALL, (cl_uint)currentPlatformDevices,
                       deviceIds.data(), nullptr);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get device id for platform \"" << platformName
                << "\". Error Code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  std::vector<cl_device_id> filteredDeviceIds;

  std::vector<std::string> filteredDeviceNames;

  std::map<std::string, size_t> countLimitMap;

  for (cl_device_id deviceId : deviceIds) {
    // filter device ids
    this->configureDevice(deviceId, devicesNode, filteredDeviceIds, filteredDeviceNames,
                          countLimitMap, useConfiguration);
  }

  if (filteredDeviceIds.size() > 0) {
    OCLPlatformWrapper platformWrapper(platformId, platformName, filteredDeviceIds,
                                       filteredDeviceNames);
    //        platforms.emplace_back(platformId, platformName,
    //        filteredDeviceIds, filteredDeviceNames);
    //        OCLPlatformWrapper &platformWrapper = *(platforms.end() - 1);
    platforms.push_back(platformWrapper);

    // create linear device list
    for (size_t deviceIndex = 0; deviceIndex < filteredDeviceIds.size(); deviceIndex++) {
      this->devices.push_back(std::make_shared<OCLDevice>(
          platformWrapper.platformId, platformWrapper.deviceIds[deviceIndex], platformName,
          platformWrapper.deviceNames[deviceIndex], platformWrapper.context,
          platformWrapper.commandQueues[deviceIndex]));
    }
  }
}

void OCLManagerMultiPlatform::configureDevice(cl_device_id deviceId, json::Node &devicesNode,
                                              std::vector<cl_device_id> &filteredDeviceIds,
                                              std::vector<std::string> &filteredDeviceNames,
                                              std::map<std::string, size_t> &countLimitMap,
                                              bool useConfiguration) {
  cl_int err;

  char deviceName[128] = {0};
  err = clGetDeviceInfo(deviceId, CL_DEVICE_NAME, 128 * sizeof(char), &deviceName, nullptr);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Failed to read the device name for device: " << deviceId
                << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  if (verbose) {
    std::cout << "OCL Info: detected device, name: \"" << deviceName << "\"" << std::endl;
  }

  if (useConfiguration) {
    if (!devicesNode.contains(deviceName)) {
      return;
    }
    if (countLimitMap.count(deviceName) > 0) {
      if (countLimitMap[deviceName] >= devicesNode[deviceName]["COUNT"].getUInt()) {
        return;
      } else {
        countLimitMap[deviceName] += 1;
      }
    } else if (devicesNode[deviceName].contains("COUNT")) {
      if (countLimitMap[deviceName] >= devicesNode[deviceName]["COUNT"].getUInt()) {
        return;
      } else {
        countLimitMap[deviceName] = 1;
      }
    }
  } else {
    if (!devicesNode.contains(deviceName)) {
      json::Node &deviceNode = devicesNode.addDictAttr(deviceName);
      deviceNode.addDictAttr("KERNELS");
    }
  }

  if (verbose) {
    std::cout << "OCL Info: using device, name: \"" << deviceName << "\"" << std::endl;
  }

  filteredDeviceIds.push_back(deviceId);
  filteredDeviceNames.push_back(deviceName);
  overallDeviceCount += 1;
}

std::vector<std::shared_ptr<OCLDevice>> &OCLManagerMultiPlatform::getDevices() {
  return this->devices;
}
}  // namespace base
}  // namespace SGPP

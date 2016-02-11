// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>

#include <sgpp/base/exception/operation_exception.hpp>

#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <vector>

namespace SGPP {
namespace base {

OCLManagerMultiPlatform::OCLManagerMultiPlatform(bool verbose) {
  parameters = std::make_shared<OCLOperationConfiguration>();
  parameters->replaceIDAttr("VERBOSE", verbose);
  parameters->replaceIDAttr("OCL_MANAGER_VERBOSE", false);
  parameters->replaceIDAttr("SHOW_BUILD_LOG", false);
  parameters->replaceDictAttr("PLATFORMS");
  parameters->replaceIDAttr("LOAD_BALANCING_VERBOSE", false);

  this->verbose = verbose;
  overallDeviceCount = 0;

  this->configure(*parameters, false);
}

OCLManagerMultiPlatform::OCLManagerMultiPlatform(
  std::shared_ptr<base::OCLOperationConfiguration> parameters) :
  parameters(parameters) {
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
        errorString << "OCL Error: Could not release command queue! "
                    "Error Code: " << err
                    << std::endl;
        throw SGPP::base::operation_exception(errorString.str());
      }
    }

    cl_int err = clReleaseContext(platform.context);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Could not release context! "
                  "Error Code: " << err <<
                  std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }
  }
}

void OCLManagerMultiPlatform::buildKernel(const std::string& program_src,
    const char* kernel_name,
    std::map<cl_platform_id, std::vector<cl_kernel> >& kernels) {
  cl_int err;

  for (OCLPlatformWrapper& platform : this->platforms) {
    // setting the program
    const char* kernel_src = program_src.c_str();
    cl_program program = clCreateProgramWithSource(
                           platform.context, 1, &kernel_src,
                           NULL, &err);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create program! "
                  "Error Code: " << err <<
                  std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }

    std::string build_opts;

    if (!(*parameters).contains("ENABLE_OPTIMIZATIONS")
        || (*parameters)["ENABLE_OPTIMIZATIONS"].getBool()) {
      // TODO(pfandedd): user should be able to change
      std::string optimizationFlags = "";

      if ((*parameters).contains("OPTIMIZATION_FLAGS")) {
        optimizationFlags = (*parameters)["OPTIMIZATION_FLAGS"].get();
      }

      build_opts = optimizationFlags;
      // -O5  -cl-mad-enable -cl-denorms-are-zero -cl-no-signed-zeros
      // -cl-unsafe-math-optimizations -cl-finite-math-only
      // -cl-fast-relaxed-math
    } else {
      build_opts = "-cl-opt-disable";  // -g
    }

    // compiling the program
    err = clBuildProgram(program, 0, NULL, build_opts.c_str(), NULL, NULL);

    if (err != CL_SUCCESS) {
      // get the build log
      size_t len;
      clGetProgramBuildInfo(program, platform.deviceIds[0],
                            CL_PROGRAM_BUILD_LOG, 0,
                            NULL, &len);
      std::string buffer(len, '\0');
      clGetProgramBuildInfo(program, platform.deviceIds[0],
                            CL_PROGRAM_BUILD_LOG, len,
                            &buffer[0], NULL);
      buffer = buffer.substr(0, buffer.find('\0'));

      if (verbose) {
        std::cout << "--- Build Log ---" << std::endl << buffer << std::endl;
      }

      std::stringstream errorString;
      errorString << "OCL Error: OpenCL build error. Error code: " << err <<
                  std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }

    for (size_t i = 0; i < platform.deviceIds.size(); i++) {
      // creating the kernel
      cl_kernel kernel = clCreateKernel(program, kernel_name, &err);

      if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to create kernel! "
                    "Error code: " << err <<
                    std::endl;
        throw SGPP::base::operation_exception(errorString.str());
      }

      kernels[platform.platformId].push_back(kernel);
    }

    if (program) {
      clReleaseProgram(program);
    }
  }
}

std::shared_ptr<base::OCLOperationConfiguration>
OCLManagerMultiPlatform::getConfiguration() {
  return this->parameters;
}

void OCLManagerMultiPlatform::configure(base::OCLOperationConfiguration&
                                        configuration, bool useConfiguration) {
  cl_int err;

  // determine number of available OpenCL platforms
  cl_uint platformCount;
  err = clGetPlatformIDs(0, nullptr, &platformCount);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString <<
                "OCL Error: Unable to get number of OpenCL platforms. "
                "Error Code: " << err <<
                std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  if (verbose) {
    std::cout << "OCL Info: " << platformCount <<
              " OpenCL platforms have been found" << std::endl;
  }

  // get available platforms
  std::vector<cl_platform_id> platformIds(platformCount);
  err = clGetPlatformIDs(platformCount, platformIds.data(), nullptr);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get Platform ID. "
                "Error Code: " << err <<
                std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  for (size_t i = 0; i < platformCount; i++) {
    this->configurePlatform(platformIds[i], *parameters, useConfiguration);
  }
}

void OCLManagerMultiPlatform::configurePlatform(cl_platform_id platformId,
    base::OCLOperationConfiguration& configuration, bool useConfiguration) {
  cl_int err;

  char platformName[128] = { 0 };
  err = clGetPlatformInfo(platformId, CL_PLATFORM_NAME, 128 * sizeof(char),
                          platformName, nullptr);

  if (CL_SUCCESS != err) {
    std::stringstream errorString;
    errorString << "OCL Error: can't get platform name!" << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  } else {
    if (verbose) {
      std::cout << "OCL Info: detected platform, name: \"" <<
                platformName << "\"" <<
                std::endl;
    }
  }

  if (verbose) {
    char vendor_name[128] = { 0 };
    err = clGetPlatformInfo(platformId, CL_PLATFORM_VENDOR, 128 * sizeof(char),
                            vendor_name, nullptr);

    if (CL_SUCCESS != err) {
      std::stringstream errorString;
      errorString << "OCL Error: Can't get platform vendor!" << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    } else {
      std::cout << "OCL Info: detected platform vendor name: \"" <<
                vendor_name <<
                "\"" << std::endl;
    }
  }

  if (useConfiguration) {
    if (!(*parameters)["PLATFORMS"].contains(platformName)) {
      return;
    }
  } else {
    // creating new configuration
    json::Node& platformNode =
      (*parameters)["PLATFORMS"].addDictAttr(platformName);
    platformNode.addDictAttr("DEVICES");
  }

  if (verbose) {
    std::cout << "OCL Info: using platform, name: \"" << platformName << "\"" <<
              std::endl;
  }

  json::Node& devicesNode = (*parameters)["PLATFORMS"][platformName]["DEVICES"];

  uint32_t currentPlatformDevices;
  // get the number of devices
  err = clGetDeviceIDs(platformId, CL_DEVICE_TYPE_ALL, 0, nullptr,
                       &currentPlatformDevices);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get device count. "
                "Error Code: " << err <<
                std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  std::vector<cl_device_id> deviceIds(currentPlatformDevices);
  err = clGetDeviceIDs(platformId, CL_DEVICE_TYPE_ALL,
                       (cl_uint) currentPlatformDevices, deviceIds.data(),
                       nullptr);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get device id for platform \"" <<
                platformName << "\". Error Code: " << err
                << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  std::vector<cl_device_id> filteredDeviceIds;

  std::map<std::string, size_t> countLimitMap;

  for (cl_device_id deviceId : deviceIds) {
    // filter device ids
    this->configureDevice(deviceId, devicesNode, filteredDeviceIds,
                          countLimitMap,
                          useConfiguration);
  }

  OCLPlatformWrapper platformWrapper(platformId, platformName,
                                     filteredDeviceIds);

  platforms.push_back(platformWrapper);
}

void OCLManagerMultiPlatform::configureDevice(cl_device_id deviceId,
    json::Node& devicesNode,
    std::vector<cl_device_id>& filteredDeviceIds,
    std::map<std::string, size_t>& countLimitMap,
    bool useConfiguration) {
  cl_int err;

  char deviceName[128] = { 0 };
  err = clGetDeviceInfo(deviceId,
                        CL_DEVICE_NAME, 128 * sizeof(char), &deviceName,
                        nullptr);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Failed to read the device name "
                "for device: " <<
                deviceId << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  if (verbose) {
    std::cout << "OCL Info: detected device, name: \"" << deviceName << "\"" <<
              std::endl;
  }

  if (useConfiguration) {
    if (!devicesNode.contains(deviceName)) {
      return;
    }

    if (countLimitMap.count(deviceName) > 0) {
      if (countLimitMap[deviceName] >=
          devicesNode[deviceName]["COUNT"].getUInt()) {
        return;
      } else {
        countLimitMap[deviceName] += 1;
      }
    } else if (devicesNode[deviceName].contains("COUNT")) {
      countLimitMap[deviceName] = 1;
    }
  } else {
    if (!devicesNode.contains(deviceName)) {
      json::Node& deviceNode = devicesNode.addDictAttr(deviceName);
      deviceNode.addDictAttr("KERNELS");
    }
  }

  if (verbose) {
    std::cout << "OCL Info: using device, name: \"" << deviceName << "\"" <<
              std::endl;
  }

  filteredDeviceIds.push_back(deviceId);
  overallDeviceCount += 1;
}
cl_kernel OCLManagerMultiPlatform::buildKernel(const std::string &source, std::shared_ptr<OCLDevice> device,
        const std::string &kernelName) {
    cl_int err;

    // setting the program
    const char* kernel_src = source.c_str();
    cl_program program = clCreateProgramWithSource(device->context, 1, &kernel_src,
    NULL, &err);

    if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to create program! Error Code: " << err << std::endl;
        throw SGPP::base::operation_exception(errorString.str());
    }

    std::string build_opts;

    if (!(*parameters).contains("ENABLE_OPTIMIZATIONS") || (*parameters)["ENABLE_OPTIMIZATIONS"].getBool()) {
        //TODO: user should be able to change
        std::string optimizationFlags = "";
        if ((*parameters).contains("OPTIMIZATION_FLAGS")) {
            optimizationFlags = (*parameters)["OPTIMIZATION_FLAGS"].get();
        }
        build_opts = optimizationFlags; // -O5  -cl-mad-enable -cl-denorms-are-zero -cl-no-signed-zeros -cl-unsafe-math-optimizations -cl-finite-math-only -cl-fast-relaxed-math
    } else {
        build_opts = "-cl-opt-disable"; // -g
    }

    // compiling the program
    err = clBuildProgram(program, 0, NULL, build_opts.c_str(), NULL, NULL);

    if (err != CL_SUCCESS) {
        // get the build log
        size_t len;
        clGetProgramBuildInfo(program, device->deviceId, CL_PROGRAM_BUILD_LOG, 0,
        NULL, &len);
        std::string buffer(len, '\0');
        clGetProgramBuildInfo(program, device->deviceId, CL_PROGRAM_BUILD_LOG, len, &buffer[0], NULL);
        buffer = buffer.substr(0, buffer.find('\0'));

        if (verbose) {
            std::cout << "--- Build Log ---" << std::endl << buffer << std::endl;
        }

        std::stringstream errorString;
        errorString << "OCL Error: OpenCL build error. Error code: " << err << std::endl;
        throw SGPP::base::operation_exception(errorString.str());

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

std::vector<std::shared_ptr<OCLDevice>> &OCLManagerMultiPlatform::getDevices() {
    return this->devices;
}


}  // namespace base
}  // namespace SGPP

/*
 OCLManagerMultiPlatform::OCLManagerMultiPlatform(std::shared_ptr<base::OCLOperationConfiguration> parameters) :
 parameters(parameters), deviceType(0) {

 // augment default values to configuration

 if (parameters->contains("VERBOSE") == false) {
 // sets the kernel to verbose
 parameters->addIDAttr("VERBOSE", true);
 }

 if (parameters->contains("OCL_MANAGER_VERBOSE") == false) {
 // sets the manager to verbose
 parameters->addIDAttr("OCL_MANAGER_VERBOSE", false);
 }

 if (parameters->contains("SHOW_BUILD_LOG") == false) {
 parameters->addIDAttr("SHOW_BUILD_LOG", false);
 }

 if (parameters->contains("PLATFORMS") == false) {
 parameters->addListAttr("PLATFORMS");
 }

 this->verbose = (*parameters)["OCL_MANAGER_VERBOSE"].getBool();

 cl_int err;

 this->setPlatformIDs();

 if (verbose) {
 this->printPlatformsInfo();
 }

 this->setupPlatforms();

 this->setDeviceType();

 this->setTotalDeviceCount();

 this->setupDeviceIDs();

 if (verbose) {
 std::cout << "OCL Info: " << overallDeviceCount << " OpenCL devices have been found!" << std::endl;
 }

 //    bool isSelectSpecificDeviceEnabled = (*parameters)["SELECT_SPECIFIC_DEVICE"].get().compare("DISABLED") != 0;
 //    bool isMaxDevicesEnabled = (*parameters)["MAX_DEVICES"].get().compare("DISABLED") != 0;
 bool isSelectSpecificDeviceEnabled = (*parameters).contains("SELECT_SPECIFIC_DEVICE");
 bool isMaxDevicesEnabled = (*parameters).contains("MAX_DEVICES");

 if (isSelectSpecificDeviceEnabled && isMaxDevicesEnabled) {
 std::stringstream errorString;
 errorString << "OCL Error: Can specify \"MAX_DEVICES\" and \"SELECT_SPECIFIC_DEVICE\" at the same time"
 << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 }

 if (isSelectSpecificDeviceEnabled) {
 if (this->platforms.size() > 1) {
 std::stringstream errorString;
 errorString << "OCL Error: Can only select a specific device if only one platform is used" << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 }

 size_t selectedDevice = (*parameters)["SELECT_SPECIFIC_DEVICE"].getUInt();

 if (selectedDevice > platforms[0].deviceCount) {
 std::stringstream errorString;
 errorString << "OCL Error: Illegal value set for \"SELECT_SPECIFIC_DEVICE\":" << selectedDevice
 << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 }

 //always for the first platform
 platforms[0].deviceIds[0] = platforms[0].deviceIds[selectedDevice];
 platforms[0].deviceCount = 1;
 overallDeviceCount = 1;

 if (verbose) {
 std::cout << "OCL Info: select device number " << selectedDevice << std::endl;
 }

 } else if (isMaxDevicesEnabled) {
 if (this->platforms.size() > 1) {
 std::stringstream errorString;
 errorString << "OCL Error: Can only select devices if only one platform is used" << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 }

 size_t maxDevices = (*parameters)["MAX_DEVICES"].getUInt();

 //always for the first platform
 platforms[0].deviceCount = maxDevices;
 overallDeviceCount = static_cast<cl_uint>(maxDevices);

 if (verbose) {
 std::cout << "OCL Info: select number of devices: " << maxDevices << std::endl;
 }
 }

 if (verbose) {
 std::cout << "OCL Info: using " << overallDeviceCount << " device/s" << std::endl;
 }

 for (OCLPlatformWrapper &platform : platforms) {
 // Create OpenCL context
 cl_context_properties properties[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties) platform.platformId, 0 };

 platform.context = clCreateContext(properties, (cl_uint) platform.getDeviceCount(), platform.deviceIds, nullptr,
 nullptr, &err);

 if (err != CL_SUCCESS) {
 std::stringstream errorString;
 errorString << "OCL Error: Failed to create OpenCL context! Error Code: " << err << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 }
 }

 // Creating the command queues
 for (OCLPlatformWrapper &platform : platforms) {
 platform.commandQueues = new cl_command_queue[platform.getDeviceCount()];
 for (size_t deviceIndex = 0; deviceIndex < platform.getDeviceCount(); deviceIndex++) {

 char buffer[128];
 err = clGetDeviceInfo(platform.deviceIds[deviceIndex],
 CL_DEVICE_NAME, 128 * sizeof(char), &buffer, nullptr);

 if (err != CL_SUCCESS) {
 std::stringstream errorString;
 errorString << "OCL Error: Failed to read the device name for device: " << deviceIndex << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 }

 if (verbose) {
 std::cout << "OCL Info: device name: " << buffer << std::endl;
 }

 platform.commandQueues[deviceIndex] = clCreateCommandQueue(platform.context,
 platform.deviceIds[deviceIndex],
 CL_QUEUE_PROFILING_ENABLE, &err);

 if (err != CL_SUCCESS) {
 std::stringstream errorString;
 errorString << "OCL Error: Failed to create command queue! Error Code: " << err << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 }
 }

 if (verbose) {
 std::cout << "OCL Info: Successfully initialized OpenCL (local workgroup size: "
 << (*parameters)["LOCAL_SIZE"].get() << ")" << std::endl << std::endl;
 }
 }
 }

 void OCLManagerMultiPlatform::setPlatformIDs() {
 cl_int err;

 // determine number of available OpenCL platforms
 cl_uint platformCount;
 err = clGetPlatformIDs(0, nullptr, &platformCount);

 if (err != CL_SUCCESS) {
 std::stringstream errorString;
 errorString << "OCL Error: Unable to get number of OpenCL platforms. Error Code: " << err << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 }

 if (verbose) {
 std::cout << "OCL Info: " << platforms.size() << " OpenCL Platforms have been found" << std::endl;
 }

 // get available platforms
 cl_platform_id *platformIds = new cl_platform_id[platformCount];
 err = clGetPlatformIDs(platformCount, platformIds, nullptr);

 for (size_t i = 0; i < platformCount; i++) {
 OCLPlatformWrapper platformWrapper(platformIds[i]);
 platforms.push_back(platformWrapper);
 }

 if (err != CL_SUCCESS) {
 std::stringstream errorString;
 errorString << "OCL Error: Unable to get Platform ID. Error Code: " << err << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 }
 delete platformIds;
 }



 void OCLManagerMultiPlatform::setupPlatforms() {
 cl_int err;

 size_t selectedPlatformIndex = 0;
 bool found = false;

 for (cl_uint i = 0; i < platforms.size(); i++) {
 char platform_name[128] = { 0 };
 err = clGetPlatformInfo(platforms[i].platformId, CL_PLATFORM_NAME, 128 * sizeof(char), platform_name, nullptr);

 if (CL_SUCCESS != err) {
 std::stringstream errorString;
 errorString << "OCL Error: Can't get platform name!" << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 } else {
 if (verbose) {
 std::cout << "OCL Info: Platform " << i << " name: " << platform_name << std::endl;
 }

 if ((*parameters)["PLATFORM"].get().compare(platform_name) == 0) {
 selectedPlatformIndex = i;
 found = true;

 if (verbose) {
 std::cout << "platform selected" << std::endl;
 }
 }
 }
 }

 if ((*parameters)["PLATFORM"].get().compare("first") == 0) {
 if (verbose) {
 std::cout << "using first platform" << std::endl;
 }
 OCLPlatformWrapper selectedPlatform = platforms[0];
 platforms.clear();
 platforms.push_back(selectedPlatform);
 } else if ((*parameters)["PLATFORM"].get().compare("all") != 0) {
 if (found) {
 OCLPlatformWrapper selectedPlatform = platforms[selectedPlatformIndex];
 platforms.clear();
 platforms.push_back(selectedPlatform);
 } else {
 std::stringstream errorString;
 errorString << "OCL Error: selected platform not found" << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 }
 }

 }

 void OCLManagerMultiPlatform::setTotalDeviceCount() {
 cl_int err;
 overallDeviceCount = 0;

 for (size_t i = 0; i < platforms.size(); i++) {
 cl_platform_id platformId = platforms[i].platformId;
 uint32_t currentPlatformDevices;
 // get the number of devices
 err = clGetDeviceIDs(platformId, this->deviceType, 0, nullptr, &currentPlatformDevices);

 if (err != CL_SUCCESS) {
 std::stringstream errorString;
 errorString << "OCL Error: Unable to get device count. Error Code: " << err << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 }

 overallDeviceCount += currentPlatformDevices;
 platforms[i].deviceCount = currentPlatformDevices;
 }
 }

 void OCLManagerMultiPlatform::setupDeviceIDs() {
 cl_int err;

 for (size_t i = 0; i < platforms.size(); i++) {
 cl_platform_id platformId = platforms[i].platformId;
 cl_device_id *deviceIds = new cl_device_id[platforms[i].deviceCount];

 err = clGetDeviceIDs(platformId, this->deviceType, (cl_uint) platforms[i].deviceCount, deviceIds, nullptr);

 platforms[i].deviceIds = deviceIds;

 if (err != CL_SUCCESS) {
 std::stringstream errorString;
 errorString << "OCL Error: Unable to get device id for platform " << i << ". Error Code: " << err
 << std::endl;
 throw SGPP::base::operation_exception(errorString.str());
 }
 }
 }

 void OCLManagerMultiPlatform::setDeviceType() {
 if ((*parameters)["DEVICE_TYPE"].get() == "CL_DEVICE_TYPE_CPU") {
 if (verbose) {
 std::cout << "OCL Info: looking for CPU device" << std::endl;
 }

 this->deviceType = CL_DEVICE_TYPE_CPU;
 } else if ((*parameters)["DEVICE_TYPE"].get() == "CL_DEVICE_TYPE_GPU") {
 if (verbose) {
 std::cout << "OCL Info: looking for GPU device" << std::endl;
 }

 this->deviceType = CL_DEVICE_TYPE_GPU;
 } else if ((*parameters)["DEVICE_TYPE"].get() == "CL_DEVICE_TYPE_ACCELERATOR") {
 if (verbose) {
 std::cout << "OCL Info: looking for device of accelerator type" << std::endl;
 }

 this->deviceType = CL_DEVICE_TYPE_ACCELERATOR;
 } else if ((*parameters)["DEVICE_TYPE"].get() == "CL_DEVICE_TYPE_ALL") {
 if (verbose) {
 std::cout << "OCL Info: looking for device of all available devices" << std::endl;
 }

 this->deviceType = CL_DEVICE_TYPE_ALL;
 } else {
 throw SGPP::base::operation_exception(
 "OCL Error: No device found or unknown type specified (supported are: CPU, GPU, accelerator and all)");
 }
 }
 */

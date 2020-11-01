// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <sstream>
#include <string>

#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

namespace sgpp {
namespace base {

OCLManager::OCLManager(std::shared_ptr<base::OCLOperationConfiguration> parameters)
    : parameters(parameters) {
  // augment default values to configuration
  if (parameters->contains("LOCAL_SIZE") == false) {
    parameters->addIDAttr("LOCAL_SIZE", UINT64_C(64));
  }

  if (parameters->contains("ENABLE_OPTIMIZATIONS") == false) {
    parameters->addIDAttr("ENABLE_OPTIMIZATIONS", true);
  }

  if (parameters->contains("OPTIMIZATION_FLAGS") == false) {
    parameters->addTextAttr("OPTIMIZATION_FLAGS", "");
  }

  if (parameters->contains("INTERNAL_PRECISION") == false) {
    parameters->addTextAttr("INTERNAL_PRECISION", "double");
  }

  if (parameters->contains("PLATFORM") == false) {
    parameters->addTextAttr("PLATFORM", "first");
  }

  if (parameters->contains("DEVICE_TYPE") == false) {
    parameters->addTextAttr("DEVICE_TYPE", "CL_DEVICE_TYPE_ALL");
  }

  if (parameters->contains("REUSE_SOURCE") == false) {
    parameters->addIDAttr("REUSE_SOURCE", false);
  }

  if (parameters->contains("WRITE_SOURCE") == false) {
    parameters->addIDAttr("WRITE_SOURCE", false);
  }

  if (parameters->contains("VERBOSE") == false) {
    // sets the kernel to verbose
    parameters->addIDAttr("VERBOSE", true);
  }

  if (parameters->contains("OCL_MANAGER_VERBOSE") == false) {
    // sets the manager to verbose
    parameters->addIDAttr("OCL_MANAGER_VERBOSE", false);
  }

  if (parameters->contains("LOAD_BALANCING_VERBOSE") == false) {
    parameters->addIDAttr("LOAD_BALANCING_VERBOSE", false);
  }

  if (parameters->contains("SHOW_BUILD_LOG") == false) {
    parameters->addIDAttr("SHOW_BUILD_LOG", false);
  }

  verbose = (*parameters)["OCL_MANAGER_VERBOSE"].getBool();
  cl_int err;

  // upper limit for number of devices of a single platform

  // determine number of available OpenCL platforms
  err = clGetPlatformIDs(0, nullptr, &num_platforms);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get number of OpenCL platforms. "
                   "Error Code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  if (verbose) {
    std::cout << "OCL Info: " << num_platforms << " OpenCL Platforms have been found" << std::endl;
  }

  // get available platforms
  platform_ids = new cl_platform_id[num_platforms];
  err = clGetPlatformIDs(num_platforms, platform_ids, nullptr);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get Platform ID. "
                   "Error Code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  for (cl_uint ui = 0; ui < num_platforms; ui++) {
    char vendor_name[128] = {0};
    err = clGetPlatformInfo(platform_ids[ui], CL_PLATFORM_VENDOR, 128 * sizeof(char), vendor_name,
                            nullptr);

    if (CL_SUCCESS != err) {
      std::stringstream errorString;
      errorString << "OCL Error: Can't get platform vendor!" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    } else {
      if (vendor_name != nullptr && verbose) {
        std::cout << "OCL Info: Platform " << ui << " vendor name: " << vendor_name << std::endl;
      }
    }

    char version_info[128] = {0};
    err = clGetPlatformInfo(platform_ids[ui], CL_PLATFORM_VERSION, 128 * sizeof(char), version_info,
                            nullptr);

    if (CL_SUCCESS != err) {
      std::stringstream errorString;
      errorString << "OCL Error: Can't get platform version!" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    } else {
      if (version_info != nullptr && verbose) {
        std::cout << "OCL Info: Platform " << ui << " version: " << version_info << std::endl;
      }
    }

    char platform_name[128] = {0};
    err = clGetPlatformInfo(platform_ids[ui], CL_PLATFORM_NAME, 128 * sizeof(char), platform_name,
                            nullptr);

    if (CL_SUCCESS != err) {
      std::stringstream errorString;
      errorString << "OCL Error: Can't get platform name!" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    } else {
      if (platform_name != nullptr && verbose) {
        std::cout << "OCL Info: Platform " << ui << " name: " << platform_name << std::endl;
      }

      if ((*parameters)["PLATFORM"].get().compare(platform_name) == 0) {
        platform_id = platform_ids[ui];

        if (verbose) {
          std::cout << "platform selected" << std::endl;
        }
      }
    }
  }

  if (verbose) {
    std::cout << std::endl;
  }

  if ((*parameters)["PLATFORM"].get().compare("first") == 0) {
    if (verbose) {
      std::cout << "using first platform" << std::endl;
    }

    platform_id = platform_ids[0];
  }

  // Find out how many devices there are
  if ((*parameters)["DEVICE_TYPE"].get() == "CL_DEVICE_TYPE_CPU") {
    if (verbose) {
      std::cout << "OCL Info: looking for CPU device" << std::endl;
    }

    // get the number of devices
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, 0, nullptr, &num_devices);
    device_ids = new cl_device_id[num_devices];
    // get the device ids
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, num_devices, device_ids, nullptr);
  } else if ((*parameters)["DEVICE_TYPE"].get() == "CL_DEVICE_TYPE_GPU") {
    if (verbose) {
      std::cout << "OCL Info: looking for GPU device" << std::endl;
    }

    // get the number of devices
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 0, nullptr, &num_devices);
    device_ids = new cl_device_id[num_devices];
    // get the device ids
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, num_devices, device_ids, nullptr);
  } else if ((*parameters)["DEVICE_TYPE"].get() == "CL_DEVICE_TYPE_ACCELERATOR") {
    if (verbose) {
      std::cout << "OCL Info: looking for device of accelerator type" << std::endl;
    }

    // get the number of devices
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ACCELERATOR, 0, nullptr, &num_devices);
    device_ids = new cl_device_id[num_devices];
    // get the device ids
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ACCELERATOR, num_devices, device_ids, nullptr);
  } else if ((*parameters)["DEVICE_TYPE"].get() == "CL_DEVICE_TYPE_ALL") {
    if (verbose) {
      std::cout << "OCL Info: looking for device of all available devices" << std::endl;
    }

    // get the number of devices
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, 0, nullptr, &num_devices);
    device_ids = new cl_device_id[num_devices];
    // get the device ids
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, num_devices, device_ids, nullptr);
  } else {
    throw sgpp::base::operation_exception(
        "OCL Error: No device found or unknown type specified "
        "(supported are: CPU, GPU, accelerator and all)");
  }

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get Device IDs. Error Code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  if (verbose) {
    std::cout << "OCL Info: " << num_devices << " OpenCL devices have been found!" << std::endl;
  }

  bool isMaxDevicesEnabled = (*parameters).contains("MAX_DEVICES");

  if ((*parameters).contains("SELECT_SPECIFIC_DEVICE")) {
    if (isMaxDevicesEnabled) {
      std::stringstream errorString;
      errorString << "OCL Error: Cannot select a specific device if more than "
                     "one device is used, MAX_DEVICES be set incorrectly" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    size_t selectedDevice = (*parameters)["SELECT_SPECIFIC_DEVICE"].getUInt();

    if (selectedDevice > num_devices) {
      std::stringstream errorString;
      errorString << "OCL Error: Illegal value set for "
                     "\"SELECT_SPECIFIC_DEVICE\"" << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    device_ids[0] = device_ids[selectedDevice];

    if (verbose) {
      std::cout << "OCL Info: select device number " << selectedDevice << std::endl;
    }
  }

  if (isMaxDevicesEnabled) {
    cl_uint maxDevices = (cl_uint)(*parameters)["MAX_DEVICES"].getUInt();

    if (maxDevices < num_devices) {
      num_devices = maxDevices;
    }
  }

  if (verbose) {
    std::cout << "OCL Info: using " << num_devices << " device/s" << std::endl;
  }

  // allocate arrays
  command_queue = new cl_command_queue[num_devices];

  // Create OpenCL context
  context = clCreateContext(0, num_devices, device_ids, nullptr, nullptr, &err);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Failed to create OpenCL context! "
                   "Error Code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  // Creating the command queues
  for (size_t i = 0; i < num_devices; i++) {
    char buffer[128];
    err = clGetDeviceInfo(device_ids[i], CL_DEVICE_NAME, 128 * sizeof(char), &buffer, nullptr);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to read the device name for "
                     "device: " << i << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }

    if (verbose) {
      std::cout << "OCL Info: device name: " << buffer << std::endl;
    }

    command_queue[i] =
        clCreateCommandQueue(context, device_ids[i], CL_QUEUE_PROFILING_ENABLE, &err);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create command queue! "
                     "Error Code: " << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  if (verbose) {
    std::cout << "OCL Info: Successfully initialized OpenCL "
                 "(local workgroup size: " << (*parameters)["LOCAL_SIZE"].getUInt() << ")"
              << std::endl
              << std::endl;
  }
}

OCLManager::~OCLManager() {
  cl_int err;

  for (size_t i = 0; i < num_devices; i++) {
    err = clReleaseCommandQueue(command_queue[i]);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Could not release command queue! "
                     "Error Code: " << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  err = clReleaseContext(context);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Could not release context! "
                   "Error Code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  delete[] this->command_queue;
}

void OCLManager::buildKernel(const std::string& program_src, const char* kernel_name,
                             cl_context context, size_t num_devices, cl_device_id* device_ids,
                             cl_kernel* kernel) {
  cl_int err;

  // setting the program
  const char* kernel_src = program_src.c_str();
  cl_program program = clCreateProgramWithSource(context, 1, &kernel_src, nullptr, &err);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Failed to create program! Error Code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  std::string build_opts;

  if ((*parameters)["ENABLE_OPTIMIZATIONS"].getBool()) {
    build_opts = (*parameters)["OPTIMIZATION_FLAGS"].get();
    // -O5  -cl-mad-enable -cl-denorms-are-zero -cl-no-signed-zeros
    // -cl-unsafe-math-optimizations -cl-finite-math-only -cl-fast-relaxed-math
  } else {
    build_opts = "-cl-opt-disable";  // -g
  }

  // compiling the program
  err = clBuildProgram(program, 0, nullptr, build_opts.c_str(), nullptr, nullptr);

  // if (err != CL_SUCCESS) {

  if ((*parameters)["SHOW_BUILD_LOG"].getBool()) {
    // get the build log
    size_t len;
    clGetProgramBuildInfo(program, device_ids[0], CL_PROGRAM_BUILD_LOG, 0, nullptr, &len);
    std::string buffer(len, '\0');
    clGetProgramBuildInfo(program, device_ids[0], CL_PROGRAM_BUILD_LOG, len, &buffer[0], nullptr);
    buffer = buffer.substr(0, buffer.find('\0'));

    if (verbose) {
      std::cout << "--- Build Log ---" << std::endl
                << buffer << std::endl;
    }
  }

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: OpenCL build error. Error code: " << err << std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  // creating the kernel
  for (size_t i = 0; i < num_devices; i++) {
    kernel[i] = clCreateKernel(program, kernel_name, &err);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel! "
                     "Error code: " << err << std::endl;
      throw sgpp::base::operation_exception(errorString.str());
    }
  }

  if (program) {
    clReleaseProgram(program);
  }
}

}  // namespace base
}  // namespace sgpp

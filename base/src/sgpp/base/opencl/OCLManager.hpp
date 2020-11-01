// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

// define required for clCreateCommandQueue on platforms
// that don't support OCL2.0 yet
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#include <CL/cl.h>
#include <string>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

namespace sgpp {
namespace base {

class OCLManager {
 public:
  std::shared_ptr<base::OCLOperationConfiguration> parameters;
  cl_uint num_platforms;
  cl_platform_id platform_id;
  cl_platform_id* platform_ids;
  cl_device_id* device_ids;
  cl_uint num_devices;
  cl_command_queue* command_queue;
  cl_context context;
  bool verbose;

 public:
  explicit OCLManager(
      std::shared_ptr<base::OCLOperationConfiguration> parameters);

  ~OCLManager();

  /**
   * @brief buildKernel builds the program that is represented by program_src
   * and creates num_devices kernel objects
   * that are stored into the array @a kernel (must be already allocated with at
   * least num_devices )
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
   */
  void buildKernel(const std::string& program_src, const char* kernel_name,
                   cl_context context, size_t num_devices,
                   cl_device_id* device_ids, cl_kernel* kernel);
};

}  // namespace base
}  // namespace sgpp

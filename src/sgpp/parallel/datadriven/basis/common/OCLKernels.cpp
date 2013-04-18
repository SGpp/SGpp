/******************************************************************************
* Copyright (C) 2010-2013 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/basis/common/OCLKernels.hpp"

#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>

// include OpenCL
#include "CL/cl.h"

namespace sg {

  namespace parallel {

    cl_int err;
    cl_platform_id platform_id;
    cl_platform_id* platform_ids;
    cl_device_id* device_ids;
    cl_uint num_platforms;
    cl_uint num_devices;
    cl_context context;
    cl_command_queue command_queue[MAX_OCL_DEVICE_COUNT];


    cl_mem clDataDP[MAX_OCL_DEVICE_COUNT], clLevelDP[MAX_OCL_DEVICE_COUNT], clIndexDP[MAX_OCL_DEVICE_COUNT];

    cl_kernel kernel_multTransDP[MAX_OCL_DEVICE_COUNT];
    cl_program program_multTransDP;

    cl_kernel kernel_multDP[MAX_OCL_DEVICE_COUNT];
    cl_program program_multDP;

    cl_kernel kernel_multTransModDP[MAX_OCL_DEVICE_COUNT];
    cl_program program_multTransModDP;

    cl_kernel kernel_multModDP[MAX_OCL_DEVICE_COUNT];
    cl_program program_multModDP;


    cl_mem clDataSP[MAX_OCL_DEVICE_COUNT], clLevelSP[MAX_OCL_DEVICE_COUNT], clIndexSP[MAX_OCL_DEVICE_COUNT];

    cl_kernel kernel_multTransSP[MAX_OCL_DEVICE_COUNT];
    cl_program program_multTransSP;

    cl_kernel kernel_multSP[MAX_OCL_DEVICE_COUNT];
    cl_program program_multSP;

    cl_kernel kernel_multTransModSP[MAX_OCL_DEVICE_COUNT];
    cl_program program_multTransModSP;

    cl_kernel kernel_multModSP[MAX_OCL_DEVICE_COUNT];
    cl_program program_multModSP;

    OCLKernels::OCLKernels() {
      // determine number of available OpenCL platforms
      err = clGetPlatformIDs(0, NULL, &num_platforms);

      if (err != CL_SUCCESS) {
        std::cout << "OCL Error: Unable to get number of OpenCL platforms. Error Code: " << err << std::endl;
      }

      std::cout << "OCL Info: " << num_platforms << " OpenCL Platforms have been found" << std::endl;

      // get available platforms
      platform_ids = new cl_platform_id[num_platforms];
      err = clGetPlatformIDs(num_platforms, platform_ids, NULL);

      if (err != CL_SUCCESS) {
        std::cout << "OCL Error: Unable to get Platform ID. Error Code: " << err << std::endl;
      }

      for (cl_uint ui = 0; ui < num_platforms; ui++) {
        char vendor_name[128] = {0};
        err = clGetPlatformInfo(platform_ids[ui], CL_PLATFORM_VENDOR, 128 * sizeof(char), vendor_name, NULL);

        if (CL_SUCCESS != err) {
          std::cout << "OCL Error: Can't get platform vendor!" << std::endl;
        } else {
          if (vendor_name != NULL) {
            std::cout << "OCL Info: Platform " << ui << " vendor name: " << vendor_name << std::endl;
          }

#ifdef USEOCL_INTEL

          if (strcmp(vendor_name, "Intel(R) Corporation") == 0) {
#ifdef USEOCL_CPU
            std::cout << "OCL Info: Using CPU Platform: " << vendor_name << std::endl;
#else
            std::cout << "OCL Info: Using GPU Platform: " << vendor_name << std::endl;
#endif
            platform_id = platform_ids[ui];
          }

#endif
#ifdef USEOCL_AMD

          if (strcmp(vendor_name, "Advanced Micro Devices, Inc.") == 0) {
#ifdef USEOCL_CPU
            std::cout << "OCL Info: Using CPU Platform: " << vendor_name << std::endl;
#else
            std::cout << "OCL Info: Using GPU Platform: " << vendor_name << std::endl;
#endif
            platform_id = platform_ids[ui];
          }

#endif
#ifdef USEOCL_NVIDIA

          if (strcmp(vendor_name, "NVIDIA Corporation") == 0) {
            std::cout << "OCL Info: Using GPU Platform: " << vendor_name << std::endl;
            platform_id = platform_ids[ui];
          }

#endif
        }
      }

      std::cout << std::endl;

      // Find out how many devices there are
#ifdef USEOCL_INTEL
      device_ids = new cl_device_id[1];
#ifdef USEOCL_CPU
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, 1, device_ids, NULL);
#else
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, device_ids, NULL);
#endif

      if (err != CL_SUCCESS) {
        std::cout << "OCL Error: Unable to get Device ID. Error Code: " << err << std::endl;
      }

      num_devices = 1;
#endif
#ifdef USEOCL_AMD
      device_ids = new cl_device_id[MAX_OCL_DEVICE_COUNT];
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, MAX_OCL_DEVICE_COUNT, device_ids, &num_devices);

      if (num_devices == 0) {
        std::cout << "OCL Error: NO GPU OpenCL devices have been found!" << std::endl;
      }

      // set max number of devices
      if (num_devices > MAX_OCL_DEVICE_COUNT) {
        num_devices = MAX_OCL_DEVICE_COUNT;
      }

#ifdef USEOCL_CPU
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, 2, device_ids, NULL);
#else
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 2, device_ids, NULL);
#endif

      if (err != CL_SUCCESS) {
        std::cout << "OCL Error: Unable to get Device ID. Error Code: " << err << std::endl;
      }

#endif
#ifdef USEOCL_NVIDIA
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, MAX_OCL_DEVICE_COUNT, NULL, &num_devices);

      if (num_devices == 0) {
        std::cout << "OCL Error: NO GPU OpenCL devices have been found!" << std::endl;
      }

      // set max number of devices
      if (num_devices > MAX_OCL_DEVICE_COUNT) {
        num_devices = MAX_OCL_DEVICE_COUNT;
      }

      device_ids = new cl_device_id[num_devices];
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, num_devices, device_ids, NULL);

      if (err != CL_SUCCESS) {
        std::cout << "OCL Error: Unable to get Device ID. Error Code: " << err << std::endl;
      }

#endif
      std::cout << "OCL Info: " << num_devices << " OpenCL devices have been found!" << std::endl;

      // Create OpenCL context
      context = clCreateContext(0, num_devices, device_ids, NULL, NULL, &err);

      if (err != CL_SUCCESS) {
        std::cout << "OCL Error: Failed to create OpenCL context! Error Code: " << err << std::endl;
      }

      // Creating the command queues
      for (size_t i = 0; i < num_devices; i++) {
#ifdef USEOCL_CPU
        command_queue[i] = clCreateCommandQueue(context, device_ids[i], CL_QUEUE_PROFILING_ENABLE , &err);
#else
        command_queue[i] = clCreateCommandQueue(context, device_ids[i], CL_QUEUE_PROFILING_ENABLE, &err);
#endif

        if (err != CL_SUCCESS) {
          std::cout << "OCL Error: Failed to create command queue! Error Code: " << err << std::endl;
        }
      }

      std::cout << std::endl;

      isFirstTimeMultTransSP = true;
      isFirstTimeMultSP = true;
      isFirstTimeMultTransModSP = true;
      isFirstTimeMultModSP = true;
      isVeryFirstTimeSP = true;

      isFirstTimeMultTransDP = true;
      isFirstTimeMultDP = true;
      isFirstTimeMultTransModDP = true;
      isFirstTimeMultModDP = true;
      isVeryFirstTimeDP = true;
    }

    OCLKernels::~OCLKernels() {
      if (!isFirstTimeMultSP) {
        for (size_t i = 0; i < num_devices; i++) {
          clReleaseKernel(kernel_multSP[i]);
        }

        clReleaseProgram(program_multSP);
      }

      if (!isFirstTimeMultDP) {
        for (size_t i = 0; i < num_devices; i++) {
          clReleaseKernel(kernel_multDP[i]);
        }

        clReleaseProgram(program_multDP);
      }

      if (!isFirstTimeMultModSP) {
        for (size_t i = 0; i < num_devices; i++) {
          clReleaseKernel(kernel_multModSP[i]);
        }

        clReleaseProgram(program_multModSP);
      }

      if (!isFirstTimeMultModDP) {
        for (size_t i = 0; i < num_devices; i++) {
          clReleaseKernel(kernel_multModDP[i]);
        }

        clReleaseProgram(program_multModDP);
      }

      if (!isFirstTimeMultTransSP) {
        for (size_t i = 0; i < num_devices; i++) {
          clReleaseMemObject(clLevelSP[i]);
          clReleaseMemObject(clIndexSP[i]);
          clReleaseKernel(kernel_multTransSP[i]);
        }

        clReleaseProgram(program_multTransSP);
      }

      if (!isFirstTimeMultTransDP) {
        for (size_t i = 0; i < num_devices; i++) {
          clReleaseMemObject(clLevelDP[i]);
          clReleaseMemObject(clIndexDP[i]);
          clReleaseKernel(kernel_multTransDP[i]);
        }

        clReleaseProgram(program_multTransDP);
      }

      if (!isFirstTimeMultTransModSP) {
        for (size_t i = 0; i < num_devices; i++) {
          clReleaseMemObject(clLevelSP[i]);
          clReleaseMemObject(clIndexSP[i]);
          clReleaseKernel(kernel_multTransModSP[i]);
        }

        clReleaseProgram(program_multTransModSP);
      }

      if (!isFirstTimeMultTransModDP) {
        for (size_t i = 0; i < num_devices; i++) {
          clReleaseMemObject(clLevelDP[i]);
          clReleaseMemObject(clIndexDP[i]);
          clReleaseKernel(kernel_multTransModDP[i]);
        }

        clReleaseProgram(program_multTransModDP);
      }

      if (!isVeryFirstTimeSP) {
        for (size_t i = 0; i < num_devices; i++) {
          clReleaseMemObject(clDataSP[i]);
        }
      }

      if (!isVeryFirstTimeDP) {
        for (size_t i = 0; i < num_devices; i++) {
          clReleaseMemObject(clDataDP[i]);
        }
      }

      for (size_t i = 0; i < num_devices; i++) {
        clReleaseCommandQueue(command_queue[i]);
      }

      clReleaseContext(context);

      delete[] device_ids;
    }

#include "parallel/datadriven/basis/common/OCLKernels_DP.imp"
#include "parallel/datadriven/basis/common/OCLKernels_SP.imp"
#include "parallel/datadriven/basis/common/OCLKernels_ModDP.imp"
#include "parallel/datadriven/basis/common/OCLKernels_ModSP.imp"

    void OCLKernels::resetKernels() {
      // linear float
      if (!isFirstTimeMultSP) {
        clReleaseProgram(program_multSP);

        for (size_t i = 0; i < num_devices; i++) {
          clReleaseKernel(kernel_multSP[i]);
        }

        isFirstTimeMultSP = true;
      }

      if (!isFirstTimeMultTransSP) {
        clReleaseProgram(program_multTransSP);

        for (size_t i = 0; i < num_devices; i++) {
          clReleaseMemObject(clLevelSP[i]);
          clReleaseMemObject(clIndexSP[i]);
          clReleaseKernel(kernel_multTransSP[i]);
        }

        isFirstTimeMultTransSP = true;
      }

      // linear double
      if (!isFirstTimeMultDP) {
        clReleaseProgram(program_multDP);

        for (size_t i = 0; i < num_devices; i++) {
          clReleaseKernel(kernel_multDP[i]);
        }

        isFirstTimeMultDP = true;
      }

      if (!isFirstTimeMultTransDP) {
        clReleaseProgram(program_multTransDP);

        for (size_t i = 0; i < num_devices; i++) {
          clReleaseMemObject(clLevelDP[i]);
          clReleaseMemObject(clIndexDP[i]);
          clReleaseKernel(kernel_multTransDP[i]);
        }

        isFirstTimeMultTransDP = true;
      }

      // modlinear double
      if (!isFirstTimeMultModDP) {
        clReleaseProgram(program_multModDP);

        for (size_t i = 0; i < num_devices; i++) {
          clReleaseKernel(kernel_multModDP[i]);
        }

        isFirstTimeMultModDP = true;
      }

      if (!isFirstTimeMultTransModDP) {
        clReleaseProgram(program_multTransModDP);

        for (size_t i = 0; i < num_devices; i++) {
          clReleaseMemObject(clLevelDP[i]);
          clReleaseMemObject(clIndexDP[i]);
          clReleaseKernel(kernel_multTransModDP[i]);
        }

        isFirstTimeMultTransModDP = true;
      }

      // modlinear float
      if (!isFirstTimeMultModSP) {
        clReleaseProgram(program_multModSP);

        for (size_t i = 0; i < num_devices; i++) {
          clReleaseKernel(kernel_multModSP[i]);
        }

        isFirstTimeMultModSP = true;
      }

      if (!isFirstTimeMultTransModSP) {
        clReleaseProgram(program_multTransModSP);

        for (size_t i = 0; i < num_devices; i++) {
          clReleaseMemObject(clLevelSP[i]);
          clReleaseMemObject(clIndexSP[i]);
          clReleaseKernel(kernel_multTransModSP[i]);
        }

        isFirstTimeMultTransModSP = true;
      }
    }

    void OCLKernels::resetData() {
      if (!isVeryFirstTimeSP) {
        for (size_t i = 0; i < num_devices; i++) {
          clReleaseMemObject(clDataSP[i]);
        }

        isVeryFirstTimeSP = true;
      }

      if (!isVeryFirstTimeDP) {
        for (size_t i = 0; i < num_devices; i++) {
          clReleaseMemObject(clDataDP[i]);
        }

        isVeryFirstTimeDP = true;
      }
    }

  }

}

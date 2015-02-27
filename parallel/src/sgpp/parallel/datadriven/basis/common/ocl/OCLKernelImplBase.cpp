// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USEOCL
#include <sgpp/parallel/datadriven/basis/common/ocl/OCLKernelImplBase.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    unsigned int get_ocl_local_size() {
      // read environment variable for local work group size: SGPP_OCL_LOCAL_SIZE
      const char* ocl_local_size_env = getenv("SGPP_OCL_LOCAL_SIZE");

      if (ocl_local_size_env != NULL) {
        unsigned int num_ocl_devices_envvalue = (unsigned int)(strtoul (ocl_local_size_env, NULL, 0));

        if (num_ocl_devices_envvalue != 0) {
          return num_ocl_devices_envvalue;
        } else {
          std::cout << "Ignoring value: \"" << ocl_local_size_env << "\" for SGPP_OCL_LOCAL_SIZE" << std::endl;
        }
      }

      return 64;
    }

    unsigned int OCLKernelImplBase::ocl_local_size = get_ocl_local_size();

    OCLKernelImplBase::OCLKernelImplBase() {
      // read number of OpenCL devices environment variable: SGPP_NUM_OCL_DEVICES
      const char* num_ocl_devices_env = getenv("SGPP_NUM_OCL_DEVICES");
      unsigned int max_number_ocl_devices = std::numeric_limits<unsigned int>::max();

      if (num_ocl_devices_env != NULL) {
        unsigned int num_ocl_devices_limit = (unsigned int)(strtoul (num_ocl_devices_env, NULL, 0));

        if (num_ocl_devices_limit != 0) {
          max_number_ocl_devices = std::min<unsigned int>(max_number_ocl_devices, num_ocl_devices_limit);
        } else {
          std::cout << "Ignoring value: \"" << num_ocl_devices_env << "\" for SGPP_NUM_OCL_DEVICES" << std::endl;
        }
      } else {
        max_number_ocl_devices = 8;
      }

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
#ifdef USEOCL_MIC
            std::cout << "OCL Info: Using MIC Platform: " << vendor_name << std::endl;
#elif USEOCL_CPU
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
#ifdef USEOCL_MIC
      device_ids = new cl_device_id[max_number_ocl_devices];
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ACCELERATOR, max_number_ocl_devices, device_ids, &num_devices);

      if (num_devices == 0) {
        std::cout << "OCL Error: NO Accelerator OpenCL devices have been found!" << std::endl;
      }

      // set max number of devices
      if (num_devices > max_number_ocl_devices) {
        num_devices = max_number_ocl_devices;
      }

#else
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
#endif
#ifdef USEOCL_AMD
      device_ids = new cl_device_id[max_number_ocl_devices];
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, max_number_ocl_devices, device_ids, &num_devices);

      if (num_devices == 0) {
        std::cout << "OCL Error: NO GPU OpenCL devices have been found!" << std::endl;
      }

      // set max number of devices
      if (num_devices > max_number_ocl_devices) {
        num_devices = max_number_ocl_devices;
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
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, max_number_ocl_devices, NULL, &num_devices);

      if (num_devices == 0) {
        std::cout << "OCL Error: NO GPU OpenCL devices have been found!" << std::endl;
      }

      // set max number of devices
      if (num_devices > max_number_ocl_devices) {
        num_devices = max_number_ocl_devices;
      }

      device_ids = new cl_device_id[num_devices];
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, num_devices, device_ids, NULL);

      if (err != CL_SUCCESS) {
        std::cout << "OCL Error: Unable to get Device ID. Error Code: " << err << std::endl;
      }

#endif
      std::cout << "OCL Info: " << num_devices << " OpenCL devices have been found!" << std::endl;

      // allocate arrays
      command_queue = new cl_command_queue[num_devices];

      clData = new cl_mem[num_devices];
      clLevel = new cl_mem[num_devices];
      clIndex = new cl_mem[num_devices];
      clMask = new cl_mem[num_devices];
      clOffset = new cl_mem[num_devices];

      clDevGrid = new cl_mem[num_devices];
      clDevTmp = new cl_mem[num_devices];
      clPinnedGrid = NULL;
      clPinnedTmp = NULL;

      kernel_mult = new cl_kernel[num_devices];
      kernel_multTrans = new cl_kernel[num_devices];

      // initialize arrays
      for (size_t i = 0; i < num_devices; i++) {
        command_queue[i] = NULL;

        clData[i] = NULL;
        clLevel[i] = NULL;
        clIndex[i] = NULL;
        clMask[i] = NULL;
        clOffset[i] = NULL;

        clDevGrid[i] = NULL;
        clDevTmp[i] = NULL;

        kernel_mult[i] = NULL;
        kernel_multTrans[i] = NULL;
      }

      // Create OpenCL context
      context = clCreateContext(0, num_devices, device_ids, NULL, NULL, &err);

      if (err != CL_SUCCESS) {
        std::cout << "OCL Error: Failed to create OpenCL context! Error Code: " << err << std::endl;
      }

      // Creating the command queues
      for (size_t i = 0; i < num_devices; i++) {
        // TODO FIXME whats the difference here?
#ifdef USEOCL_CPU
        command_queue[i] = clCreateCommandQueue(context, device_ids[i], CL_QUEUE_PROFILING_ENABLE, &err);
#else
        command_queue[i] = clCreateCommandQueue(context, device_ids[i], CL_QUEUE_PROFILING_ENABLE, &err);
#endif

        if (err != CL_SUCCESS) {
          std::cout << "OCL Error: Failed to create command queue! Error Code: " << err << std::endl;
        }
      }

      std::cout << "OCL Info: Successfully initialized OpenCL (local workgroup size: " << ocl_local_size << ")" << std::endl << std::endl;
    }

    OCLKernelImplBase::~OCLKernelImplBase() {
      releaseDataBuffers();
      releaseGridBuffers();
      releaseKernelsAndPrograms();

      // release command queue
      for (size_t i = 0; i < num_devices; i++) {
        if (command_queue[i]) {
          clReleaseCommandQueue(command_queue[i]);
        }
      }

      // release context
      clReleaseContext(context);

      delete[] command_queue;

      delete[] clData;
      delete[] clLevel;
      delete[] clIndex;
      delete[] clMask;
      delete[] clOffset;

      delete[] clDevGrid;
      delete[] clDevTmp;

      delete[] kernel_mult;
      delete[] kernel_multTrans;

      delete[] device_ids;
    }

    void OCLKernelImplBase::releaseGridBuffers() {
      for (size_t i = 0; i < num_devices; i++) {
        if (clLevel[i]) {
          clReleaseMemObject(clLevel[i]);
          clLevel[i] = NULL;
        }

        if (clIndex[i]) {
          clReleaseMemObject(clIndex[i]);
          clIndex[i] = NULL;
        }

        if (clMask[i]) {
          clReleaseMemObject(clMask[i]);
          clMask[i] = NULL;
        }

        if (clOffset[i]) {
          clReleaseMemObject(clOffset[i]);
          clOffset[i] = NULL;
        }

        if (clDevGrid[i]) {
          clReleaseMemObject(clDevGrid[i]);
          clDevGrid[i] = NULL;
        }
      }

      if (clPinnedGrid) {
        clReleaseMemObject(clPinnedGrid);
        clPinnedGrid = NULL;
      }
    }

    void OCLKernelImplBase::releaseDataBuffers() {
      for (size_t i = 0; i < num_devices; i++) {
        if (clData[i]) {
          clReleaseMemObject(clData[i]);
          clData[i] = NULL;
        }

        if (clDevTmp[i]) {
          clReleaseMemObject(clDevTmp[i]);
          clDevTmp[i] = NULL;
        }
      }

      if (clPinnedTmp) {
        clReleaseMemObject(clPinnedTmp);
        clPinnedTmp = NULL;
      }
    }

    void OCLKernelImplBase::releaseKernelsAndPrograms() {
      for (size_t i = 0; i < num_devices; i++) {
        if (kernel_mult[i]) {
          clReleaseKernel(kernel_mult[i]);
          kernel_mult[i] = NULL;
        }

        if (kernel_multTrans[i]) {
          clReleaseKernel(kernel_multTrans[i]);
          kernel_multTrans[i] = NULL;
        }
      }
    }

    void OCLKernelImplBase::resetKernel() {
      releaseGridBuffers();
    }
  }
}

#endif
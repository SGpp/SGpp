#include "OCLPDEKernels.hpp"


namespace sg {
  namespace parallel {
    namespace oclpdekernels {

      cl_uint num_devices;
      cl_uint num_platforms;
      cl_platform_id platform_id;
      cl_platform_id* platform_ids;
      cl_device_id* device_ids;
      cl_context context;
      cl_command_queue command_queue[NUMDEVS];
      size_t isVeryFirstTime = 1;
      size_t isCleanedUp = 0;
      double AverageGFLOPS = 0.0;
      // Theoretical maximum
      double GPUNVDMAXFLOPS = 660 * 0.75;
      double GPUNVDMAXFLOPSHALF = 660 * 0.375;
      size_t padding_size;
      size_t dims;
      REAL TimestepCoeff = 0.0;
      size_t LastRunBound = 1;

      void PrintGFLOPS();



      void oclCheckErr(cl_int err, const char* function) {
        if (err != CL_SUCCESS) {
          printf("Error: Failure %s: %d\n", function, err);
          exit(-1);
        }
      }



      void StartUpGPU() {
#ifdef USE_MPI
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

#ifdef USE_MPI

        if (myrank == 0) {
#endif
          std::cout << "StartUpDevice" << std::endl;
#ifdef USE_MPI
        }

#endif
        // read number of OpenCL devices environment variable: SGPP_NUM_OCL_DEVICES
        const char* num_ocl_devices_env = getenv("SGPP_NUM_OCL_DEVICES");
        unsigned int max_number_ocl_devices = std::min<unsigned int>(std::numeric_limits<unsigned int>::max(), NUMDEVS);

        if (num_ocl_devices_env != NULL) {
          unsigned int num_ocl_devices_limit = (unsigned int)(strtoul (num_ocl_devices_env, NULL, 0));

          if (num_ocl_devices_limit != 0) {
            max_number_ocl_devices = std::min<unsigned int>(max_number_ocl_devices, num_ocl_devices_limit);
          } else {
            std::cout << "Ignoring value: \"" << num_ocl_devices_env << "\" for SGPP_NUM_OCL_DEVICES" << std::endl;
          }
        } else {
          max_number_ocl_devices = std::min<unsigned int>(max_number_ocl_devices, 8);
        }

        cl_int err = CL_SUCCESS;
        err |= clGetPlatformIDs(0, NULL, &num_platforms);
        oclCheckErr(err, "OCL Error: Unable to get number of OpenCL platforms. Error Code");
#ifdef USE_MPI

        if (myrank == 0) {
#endif
          std::cout << "OCL Info: " << num_platforms << " OpenCL Platforms have been initialized" << std::endl;
#ifdef USE_MPI
        }

#endif
        platform_ids = new cl_platform_id[num_platforms];
        // get available platforms

        err |= clGetPlatformIDs(num_platforms, platform_ids, NULL);
        oclCheckErr(err, "OCL Error: Unable to get Platform ID. Error Code");

        for (cl_uint ui = 0; ui < num_platforms; ui++) {
          char vendor_name[128] = {0};
          err = clGetPlatformInfo(platform_ids[ui], CL_PLATFORM_VENDOR, 128 * sizeof(char), vendor_name, NULL);
          oclCheckErr(err, "OCL Error: Can't get platform vendor!");
#ifdef USE_MPI

          if (myrank == 0) {
#endif
            std::cout << "OCL Info: Platform " << ui << " vendor name: " << vendor_name << std::endl;
#ifdef USE_MPI
          }

#endif
#ifdef USEOCL_INTEL

          if (strcmp(vendor_name, "Intel(R) Corporation") == 0) {
#ifdef USE_MPI

            if (myrank == 0) {
#endif
#ifdef USEOCL_MIC
              std::cout << "OCL Info: Using MIC Platform: " << vendor_name << std::endl;
#elif USEOCL_CPU
              std::cout << "OCL Info: Using CPU Platform: " << vendor_name << std::endl;
#else
              std::cout << "OCL Info: Using GPU Platform: " << vendor_name << std::endl;
#endif
#ifdef USE_MPI
            }

#endif
            platform_id = platform_ids[ui];
          }

#endif
#ifdef USEOCL_AMD

          if (strcmp(vendor_name, "Advanced Micro Devices, Inc.") == 0) {
#ifdef USE_MPI

            if (myrank == 0) {
#endif
#ifdef USEOCL_CPU
              std::cout << "OCL Info: Using CPU Platform: " << vendor_name << std::endl;
#else
              std::cout << "OCL Info: Using GPU Platform: " << vendor_name << std::endl;
#endif
#ifdef USE_MPI
            }

#endif
            platform_id = platform_ids[ui];
          }

#endif
#ifdef USEOCL_NVIDIA

          if (strcmp(vendor_name, "NVIDIA Corporation") == 0) {
#ifdef USE_MPI

            if (myrank == 0) {
#endif
              std::cout << "OCL Info: Using GPU Platform: " << vendor_name << std::endl;
#ifdef USE_MPI
            }

#endif
            platform_id = platform_ids[ui];
          }

#endif
        }

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
        oclCheckErr(err, "OCL Error: Unable to get Device ID. Error Code");
        num_devices = 1;
#endif
#endif

#ifdef USEOCL_AMD
        //  //platform_id = platform_ids[0];
        //  err |= clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
        //  oclCheckErr(err, "clGetDeviceIDs AMD");
        device_ids = new cl_device_id[max_number_ocl_devices];
        err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);


        if (num_devices == 0) {
          std::cout << "OCL Error: NO GPU OpenCL devices have been found!" << std::endl;
        }

        // set max number of devices
        if (num_devices > max_number_ocl_devices) {
          num_devices = max_number_ocl_devices;
        }

#ifdef USEOCL_CPU
        err |= clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, num_devices, device_ids, NULL);
#else
        err |= clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, num_devices, device_ids, NULL);
#endif
        oclCheckErr(err, "OCL Error: Unable to get Device ID. Error Code");
#endif

#ifdef USEOCL_NVIDIA
        err |= clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, max_number_ocl_devices, NULL, &num_devices);
        oclCheckErr(err, "clGetDeviceIDs1");

        if (num_devices == 0) {
          std::cout << "OCL Error: NO GPU OpenCL devices have been found!" << std::endl;
        }

        // set max number of devices
        if (num_devices > max_number_ocl_devices) {
          num_devices = max_number_ocl_devices;
        }

        device_ids = new cl_device_id[num_devices];

        err |= clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, num_devices, device_ids, NULL);
        oclCheckErr(err, "OCL Error: Unable to get Device ID. Error Code");
#endif
#ifdef USE_MPI

        if (myrank == 0) {
#endif
          std::cout << "OCL Info: " << num_devices << " OpenCL devices have been initialized!" << std::endl;
#ifdef USE_MPI
        }

#endif

        context = clCreateContext(0, num_devices, device_ids, NULL, NULL, &err);
        oclCheckErr(err, "clCreateContext");

        for (size_t i = 0; i < num_devices; i++) {
#if QUEUEPROF
          command_queue[i] = clCreateCommandQueue(context, device_ids[i],
                                                  CL_QUEUE_PROFILING_ENABLE, &err);
#else
          command_queue[i] = clCreateCommandQueue(context, device_ids[i],
                                                  0, &err);

#endif
          oclCheckErr(err, "clCreateCommandQueue");
        }

      }



      void PrintGFLOPS() {
        Timing LaplaceInnerTiming = PrintGFLOPSLaplaceInner();
        Timing LaplaceBoundTiming = PrintGFLOPSLaplaceBound();
        Timing LTwoDotInnerTiming = PrintGFLOPSLTwoDotInner();
        Timing LTwoDotBoundTiming = PrintGFLOPSLTwoDotBound();
        Timing LTwoDotLaplaceInnerTiming = PrintGFLOPSLTwoDotLaplaceInner();
        Timing LTwoDotLaplaceBoundTiming = PrintGFLOPSLTwoDotLaplaceBound();


        double TotalCalculationTime = LaplaceInnerTiming.time +
                                      LaplaceBoundTiming.time + LTwoDotInnerTiming.time +
                                      LTwoDotBoundTiming.time + LTwoDotLaplaceInnerTiming.time +
                                      LTwoDotLaplaceBoundTiming.time;
        double LaplaceInnerWeight = LaplaceInnerTiming.time / TotalCalculationTime;
        double LaplaceBoundWeight = LaplaceBoundTiming.time / TotalCalculationTime;
        double LTwoDotInnerWeight = LTwoDotInnerTiming.time / TotalCalculationTime;
        double LTwoDotBoundWeight = LTwoDotBoundTiming.time / TotalCalculationTime;
        double LTwoDotLaplaceInnerWeight = LTwoDotLaplaceInnerTiming.time / TotalCalculationTime;
        double LTwoDotLaplaceBoundWeight = LTwoDotLaplaceBoundTiming.time / TotalCalculationTime;
        AverageGFLOPS += LaplaceInnerTiming.GFLOPS * LaplaceInnerWeight;
        AverageGFLOPS += LaplaceBoundTiming.GFLOPS * LaplaceBoundWeight;
        AverageGFLOPS += LTwoDotInnerTiming.GFLOPS * LTwoDotInnerWeight;
        AverageGFLOPS += LTwoDotBoundTiming.GFLOPS * LTwoDotBoundWeight;
        AverageGFLOPS += LTwoDotLaplaceInnerTiming.GFLOPS * LTwoDotLaplaceInnerWeight;
        AverageGFLOPS += LTwoDotLaplaceBoundTiming.GFLOPS * LTwoDotLaplaceBoundWeight;

        double AverageGOPS = LaplaceInnerTiming.GOPS * LaplaceInnerWeight;
        AverageGOPS += LaplaceBoundTiming.GOPS * LaplaceBoundWeight;
        AverageGOPS += LTwoDotInnerTiming.GOPS * LTwoDotInnerWeight;
        AverageGOPS += LTwoDotBoundTiming.GOPS * LTwoDotBoundWeight;
        AverageGOPS += LTwoDotLaplaceInnerTiming.GOPS * LTwoDotLaplaceInnerWeight;
        AverageGOPS += LTwoDotLaplaceBoundTiming.GOPS * LTwoDotLaplaceBoundWeight;

        std::cout << "Total Bound+Inner time : " << TotalCalculationTime << std::endl;
        std::cout << "Average GOPS: " << AverageGOPS << std::endl;
        std::cout << "Average GFLOPS: " << AverageGFLOPS << std::endl;

      }

      double AccumulateTiming(cl_event* GPUExecution, size_t gpuid) {
        cl_int ciErrNum = CL_SUCCESS;

        cl_ulong startTime, endTime;
        ciErrNum = clGetEventProfilingInfo(GPUExecution[gpuid], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, NULL);
        oclCheckErr(ciErrNum, "clGetEventProfilingInfo1");
        ciErrNum = clGetEventProfilingInfo(GPUExecution[gpuid], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, NULL);
        oclCheckErr(ciErrNum, "clGetEventProfilingInfo2");
        return (double)(endTime - startTime);
      }

      double AccumulateWaiting(cl_event* GPUExecution, size_t gpuid) {
        cl_int ciErrNum = CL_SUCCESS;

        cl_ulong startTime, endTime;
        ciErrNum = clGetEventProfilingInfo(GPUExecution[gpuid], CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &startTime, NULL);
        oclCheckErr(ciErrNum, "clGetEventProfilingInfo3");
        ciErrNum = clGetEventProfilingInfo(GPUExecution[gpuid], CL_PROFILING_COMMAND_SUBMIT, sizeof(cl_ulong), &endTime, NULL);
        oclCheckErr(ciErrNum, "clGetEventProfilingInfo4");
        return (double)(endTime - startTime);

      }

    }
    using namespace oclpdekernels;
    void OCLPDEKernels::CleanUpGPU() {
      if (isCleanedUp == 0) {
        delete device_ids;
        delete platform_ids;

        for (unsigned int i = 0; i < num_devices; i++)  {
          clReleaseCommandQueue(command_queue[i]);
        }

        clReleaseContext(context);
        isCleanedUp = 1;
        //  CleanUpInner();
        CleanUpBound();
#if TOTALTIMING
#ifdef USE_MPI
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        if (myrank == 0) {
#endif
          PrintGFLOPS();
#ifdef USE_MPI
        }

#endif
#endif
      }
    }
  }
}

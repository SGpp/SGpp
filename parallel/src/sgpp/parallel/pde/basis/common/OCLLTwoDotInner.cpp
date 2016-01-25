// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OCLLTwoDotInner.hpp"
using namespace SGPP::base;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    namespace oclpdekernels {


      cl_kernel LTwoDotInnerKernel[NUMDEVS];
      double MultTimeLTwoDotInner = 0.0;
      double ReduTimeLTwoDotInner = 0.0;
      double CounterLTwoDotInner = 0.0;
      double LTwoDotInnerStartupTime = 0.0;
      double LTwoDotInnerExecTime = 0.0;
      double LTwoDotInnerAllReduceTime = 0.0;

      void CompileLTwoDotInner(int id, std::string kernel_src, cl_kernel* kernel) {
        cl_int err = CL_SUCCESS;
        const char* source2 = "";
        std::string l2dotfunction = InnerLTwoDotFunction();

        std::stringstream stream_program_src;
#ifdef USEOCL_NVIDIA
        stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
#endif
        stream_program_src << "#define REAL double" << std::endl;

        stream_program_src << "#define LSIZE " << LSIZE << std::endl;
        stream_program_src << l2dotfunction << std::endl;
        stream_program_src << "__kernel void multKernel(__global  REAL* ptrLevel, " << std::endl;
        stream_program_src << "			 __constant  REAL* ptrLevelIndexLevelintcon," << std::endl;
        stream_program_src << "			 __global  REAL* ptrIndex," << std::endl;
        stream_program_src << "			 __global  REAL* ptrLevel_int," << std::endl;
        stream_program_src << "			 __global  REAL* ptrAlpha," << std::endl;
        stream_program_src << "			 __global  REAL* ptrParResult," << std::endl;
        stream_program_src << "			 __constant  REAL* ptrLcl_q," << std::endl;
        stream_program_src << "			 ulong overallMultOffset," << std::endl;
        stream_program_src << "			 ulong ConstantMemoryOffset)" << std::endl;
        stream_program_src << "{" << std::endl;

        cl_uint sizec;
        clGetDeviceInfo(device_ids[0], CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, sizeof(cl_uint), &sizec, 0);

        size_t keplerUseRegister = (sizec >= 3);

        // size_t fermiUseRegister = (dims <= 5);

#ifdef USEOCL_INTEL
        keplerUseRegister = 1;
        // fermiUseRegister = 0;
#endif

        stream_program_src << "__local REAL alphaTemp[" << LSIZE << "];" << std::endl;

        stream_program_src << "ptrLevel += get_local_id(0);" << std::endl;
        stream_program_src << "ptrIndex += get_local_id(0);" << std::endl;
        stream_program_src << "ptrLevel_int += get_local_id(0);" << std::endl;

        stream_program_src << "alphaTemp[get_local_id(0)]   = ptrAlpha[get_global_id(0) + overallMultOffset +  ConstantMemoryOffset*LSIZE];" << std::endl;
        stream_program_src << " barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;

        stream_program_src << "REAL res = 0.0;" << std::endl;
        const char* typeREAL = "";

        if (keplerUseRegister) {
          const char* typeREAL2 = "REAL ";

          for (size_t d_inner = 0; d_inner < dims; d_inner++) {
            stream_program_src << typeREAL2 << "i_level_" << d_inner << " =         ptrLevel[" << d_inner* storageSizePadded << " + (get_global_id(1))*LSIZE];" << std::endl;
            stream_program_src << typeREAL2 << "i_index_" << d_inner << " =         ptrIndex[" << d_inner* storageSizePadded << " + (get_global_id(1))*LSIZE];" << std::endl;
            stream_program_src << typeREAL2 << "i_level_int_" << d_inner << " =     ptrLevel_int[" << d_inner* storageSizePadded << " + (get_global_id(1))*LSIZE];" << std::endl;
          }
        }


        stream_program_src << "for (unsigned k = 0; k < " << LSIZE << "; k++) {" << std::endl;
        typeREAL = "REAL ";
        stream_program_src << typeREAL << "element = alphaTemp[" << "k" << "];" << std::endl;

        for (size_t d_inner = 0; d_inner < dims; d_inner++) {
          if (d_inner == 0) {
            typeREAL = "REAL ";
          } else {
            typeREAL = "";
          }

          std::stringstream registerId;
          registerId << "_" << d_inner;
          std::string registerIdStr = registerId.str();

          if (!keplerUseRegister) {
            const char* typeREAL2 = "REAL ";
            stream_program_src << typeREAL2 << "i_level" << registerIdStr << " =         ptrLevel[" << d_inner* storageSizePadded << " + (get_global_id(1))*LSIZE];" << std::endl;
            stream_program_src << typeREAL2 << "i_index" << registerIdStr << " =         ptrIndex[" << d_inner* storageSizePadded << " + (get_global_id(1))*LSIZE];" << std::endl;
            stream_program_src << typeREAL2 << "i_level_int" << registerIdStr << " =     ptrLevel_int[" << d_inner* storageSizePadded << " + (get_global_id(1))*LSIZE];" << std::endl;
          }

          stream_program_src << typeREAL << "j_levelreg = ptrLevelIndexLevelintcon[(get_group_id(0)*LSIZE+k)*" << dims * 3 << " + " << d_inner * 3 << "]; " << std::endl;
          stream_program_src << typeREAL << "j_indexreg = ptrLevelIndexLevelintcon[(get_group_id(0)*LSIZE+k)*" << dims * 3 << " + " << d_inner * 3 + 1 << "]; " << std::endl;
          stream_program_src << typeREAL << "j_level_intreg = ptrLevelIndexLevelintcon[(get_group_id(0)*LSIZE+k)*" << dims * 3 << " + " << d_inner * 3 + 2 << "]; " << std::endl;

          stream_program_src << "element *= (l2dot(i_level" << registerIdStr << ",\
	                        i_index" << registerIdStr << ",\
	                        i_level_int" << registerIdStr << ",\
	                        j_levelreg,j_indexreg,j_level_intreg," << "ptrLcl_q[" << d_inner << "]));" << std::endl;

        }

        stream_program_src << "res +=  element;" << std::endl;

        stream_program_src << "}" << std::endl;
        stream_program_src << "ptrParResult[((get_group_id(0)+ConstantMemoryOffset)*get_global_size(1) + get_global_id(1))*" << LSIZE << " + get_local_id(0)] = res; " << std::endl;
        stream_program_src << "}" << std::endl;

        std::string program_src = stream_program_src.str();
        source2 = program_src.c_str();
#if PRINTOCL
        std::cout <<  source2 << std::endl;
#endif
        cl_program program = clCreateProgramWithSource(context, 1, (const char**)&source2, NULL, &err);
        oclCheckErr(err, "clCreateProgramWithSource");
        char buildOptions[256];
        int ierr = snprintf(buildOptions, sizeof(buildOptions), "-DSTORAGE=%lu -DSTORAGEPAD=%lu -DNUMGROUPS=%lu -DDIMS=%lu -cl-finite-math-only -cl-strict-aliasing -cl-fast-relaxed-math ", storageSize, storageSizePadded, num_groups, dims);

        if (ierr < 0) {
          printf("Error in Build Options");
          exit(-1);
        }

        err = clBuildProgram(program, 0, NULL, buildOptions, NULL, NULL);

        if (err != CL_SUCCESS) {
          std::cout << "OCL Error: compileLTwoDotInner OpenCL Build Error. Error Code: " << err << std::endl;

          size_t len;
          char buffer[10000];

          // get the build log
          clGetProgramBuildInfo(program, device_ids[0], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);

          std::cout << "--- Build Log ---" << std::endl << buffer << std::endl;
        }

        oclCheckErr(err, "clBuildProgram");


        kernel[id] = clCreateKernel(program, kernel_src.c_str(), &err);
        oclCheckErr(err, "clCreateKernel");

        err |= clReleaseProgram(program);
        oclCheckErr(err, "clReleaseProgram");

      } // compileLTwoDotInner

      void CompileLTwoDotInnerKernels() {
        for (unsigned int i = 0; i < num_devices; ++i) {
          CompileLTwoDotInner(i, "multKernel", LTwoDotInnerKernel);

          CompileReduceInner(i, "ReduceInnerKernel", ReduceInnerKernel);
        }
      }

      void SetArgumentsLTwoDotInner() {
        cl_int ciErrNum = CL_SUCCESS;
        int counter = 0;

        for (unsigned int i = 0; i < num_devices; ++i) {
          ciErrNum |= clSetKernelArg(LTwoDotInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLevelInner[i]);
          ciErrNum |= clSetKernelArg(LTwoDotInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLevelIndexLevelintconInner[i]);

          ciErrNum |= clSetKernelArg(LTwoDotInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrIndexInner[i]);
          ciErrNum |= clSetKernelArg(LTwoDotInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLevel_intInner[i]);

          ciErrNum |= clSetKernelArg(LTwoDotInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrAlphaInner[i]);
          ciErrNum |= clSetKernelArg(LTwoDotInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrParResultInner[i]);
          ciErrNum |= clSetKernelArg(LTwoDotInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLcl_qInner[i]);
          oclCheckErr(ciErrNum, "clSetKernelArg1 Kernel Construct");

          counter = 0;
          ciErrNum |= clSetKernelArg(ReduceInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrResultInner[i]);
          ciErrNum |= clSetKernelArg(ReduceInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrParResultInner[i]);
          counter = 0;
          oclCheckErr(ciErrNum, "clSetKernelArg1 Kernel Construct");

        }
      }

      void CleanUpLTwoDotInner() {
        if (!isFirstTimeLTwoDotInner) {
          cl_int ciErrNum = CL_SUCCESS;

          for (unsigned int i = 0; i < num_devices; i++)  {
            ciErrNum |= clReleaseKernel( LTwoDotInnerKernel[i] );
            oclCheckErr(ciErrNum, "clReleaseKernel");
          }
        }
      }

      Timing PrintGFLOPSLTwoDotInner() {
        double StorD = (double)storageSize;
        double dimsD = (double)dims;
        double LTwoDotInnerTime = (MultTimeLTwoDotInner + ReduTimeLTwoDotInner) * (1e-9);
        LTwoDotInnerTime = LTwoDotInnerExecTime;

        if (!isFirstTimeLTwoDotInner) {
          std::cout << "Overall Runtime of AllReduceLTwoDotInner: " << LTwoDotInnerAllReduceTime << std::endl;
          std::cout << "Time for GPU-calculation of LTwoDotInner: " << LTwoDotInnerTime << std::endl;
        }

        // LTwoDot Inner
        double numberOfFlop = StorD * StorD * (dimsD * (22.0 + 1) + 1);


        if (!isFirstTimeLTwoDotInner) {
          double LTwoDotInnerGFLOPS = ((numberOfFlop * CounterLTwoDotInner) / LTwoDotInnerTime) * (1e-9);
          double LTwoDotInnerOPS = StorD * StorD * (dimsD * (12.0));
          double LTwoDotInnerGOPS = ((LTwoDotInnerOPS * CounterLTwoDotInner) / LTwoDotInnerTime) * (1e-9);
          return (Timing(LTwoDotInnerGFLOPS, LTwoDotInnerGOPS, LTwoDotInnerTime));
        }

        return (Timing());
      } // PrintGFLOPSLTwoDotInner


    } // namespace parallel
  }
}
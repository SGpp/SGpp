// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OCLLTwoDotBound.hpp"
using namespace SGPP::base;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    namespace oclpdekernels {

      cl_kernel LTwoDotBoundKernel[NUMDEVS];
      double MultTimeLTwoDotBound = 0.0;
      double ReduTimeLTwoDotBound = 0.0;
      double CounterLTwoDotBound = 0.0;
      double LTwoDotBoundStartupTime = 0.0;
      double LTwoDotBoundExecTime = 0.0;
      double LTwoDotBoundAllReduceTime = 0.0;

      void CompileLTwoDotBound(int id, std::string kernel_src, cl_kernel* kernel) {
        cl_int err = CL_SUCCESS;
        const char* source2 = "";
        std::string l2dotfunction = BoundLTwoDotFunction();

        std::stringstream stream_program_src;

        stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
        stream_program_src << "#define REAL double" << std::endl;
        stream_program_src << "#define LSIZE " << LSIZE << std::endl;
        stream_program_src << l2dotfunction << std::endl;

        stream_program_src << "__kernel void multKernel(__global  REAL* ptrLevel," << std::endl;
        stream_program_src << "			 __constant  REAL* ptrLevelIndexLevelintcon," << std::endl;
        stream_program_src << "			 __global  REAL* ptrIndex, " << std::endl;
        stream_program_src << "			 __global  REAL* ptrLevel_int," << std::endl;
        stream_program_src << "			 __global  REAL* ptrAlpha," << std::endl;
        stream_program_src << "			 __global  REAL* ptrParResult," << std::endl;
        stream_program_src << "			 __constant  REAL* ptrLcl_q," << std::endl;
        stream_program_src << "			 ulong overallMultOffset," << std::endl;
        stream_program_src << "			 ulong ConstantMemoryOffset)" << std::endl;
        stream_program_src << "{" << std::endl;


        stream_program_src << "  __local REAL alphaTemp[LSIZE];" << std::endl;

        stream_program_src << "alphaTemp[get_local_id(0)]   = ptrAlpha[get_global_id(0) + overallMultOffset + ConstantMemoryOffset*LSIZE];" << std::endl;
        stream_program_src << "ptrLevel += get_local_id(0);" << std::endl;
        stream_program_src << "ptrIndex += get_local_id(0);" << std::endl;
        stream_program_src << "ptrLevel_int += get_local_id(0);" << std::endl;
        stream_program_src << "REAL res = 0.0;" << std::endl;

        cl_uint sizec;
        clGetDeviceInfo(device_ids[0], CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, sizeof(cl_uint), &sizec, 0);

        size_t keplerUseRegister = (sizec >= 3);
        // size_t fermiUseRegister = (dims <= 5);


#ifdef USEOCL_INTEL
        keplerUseRegister = 0;
        // fermiUseRegister = 1;
#endif




        if (keplerUseRegister) {
          const char* typeREAL2 = "REAL ";

          for (size_t d_inner = 0; d_inner < dims; d_inner++) {
            stream_program_src << typeREAL2 << "i_level_" << d_inner << " =         ptrLevel[" << d_inner* storageInnerSizePaddedBound << "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
            stream_program_src << typeREAL2 << "i_index_" << d_inner << " =         ptrIndex[" << d_inner* storageInnerSizePaddedBound << "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
            stream_program_src << typeREAL2 << "i_level_int_" << d_inner << " =     ptrLevel_int[" << d_inner* storageInnerSizePaddedBound << "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;

          }
        }

        stream_program_src << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
        const char* typeREAL = "";

        stream_program_src << "for (unsigned k = 0; k < " << LSIZE << "; k++) {" << std::endl;
        stream_program_src << "REAL element = alphaTemp[" << "k" << "];" << std::endl;

        for (size_t d_inner = 0; d_inner < dims; d_inner++) {
          if (d_inner == 0) {
            typeREAL = "REAL ";
          } else {
            typeREAL = "";
          }

          const char* typeREAL2 = "REAL ";

          if (!keplerUseRegister) {
            stream_program_src << typeREAL2 << "i_level_" << d_inner << " =         ptrLevel[" << d_inner* storageInnerSizePaddedBound << "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
            stream_program_src << typeREAL2 << "i_index_" << d_inner << " =         ptrIndex[" << d_inner* storageInnerSizePaddedBound << "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
            stream_program_src << typeREAL2 << "i_level_int_" << d_inner << " =     ptrLevel_int[" << d_inner* storageInnerSizePaddedBound << "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
          }

          stream_program_src << typeREAL << "j_levelreg = ptrLevelIndexLevelintcon[(get_group_id(0)*LSIZE+k)*" << dims * 3 << " + " << d_inner * 3 << "]; " << std::endl;
          stream_program_src << typeREAL << "j_indexreg = ptrLevelIndexLevelintcon[(get_group_id(0)*LSIZE+k)*" << dims * 3 << " + " << d_inner * 3 + 1 << "]; " << std::endl;
          stream_program_src << typeREAL << "j_level_intreg = ptrLevelIndexLevelintcon[(get_group_id(0)*LSIZE+k)*" << dims * 3 << " + " << d_inner * 3 + 2 << "]; " << std::endl;
          stream_program_src << "element *= (l2dot(\
             i_level_" << d_inner << ",\
             i_index_" << d_inner << ",\
             i_level_int_" << d_inner << ", \
             j_levelreg,j_indexreg,j_level_intreg,ptrLcl_q[" << d_inner << "]));" << std::endl;
        }

        stream_program_src << "res += element;" << std::endl;
        stream_program_src << "}" << std::endl;
        stream_program_src << "ptrParResult[((get_group_id(0) + ConstantMemoryOffset)*get_global_size(1) + get_global_id(1))*" << LSIZE << " + get_local_id(0)] = res; " << std::endl;
        stream_program_src << "}" << std::endl;

        std::string program_src = stream_program_src.str();
        source2 = program_src.c_str();
#if PRINTOCL
        std::cout <<  source2 << std::endl;
#endif
        cl_program program = clCreateProgramWithSource(context, 1, (const char**)&source2, NULL, &err);
        oclCheckErr(err, "clCreateProgramWithSource");
        char buildOptions[256];
        int ierr = snprintf(buildOptions, sizeof(buildOptions), "-DSTORAGE=%lu -DSTORAGEPAD=%lu -DNUMGROUPS=%lu -DDIMS=%lu -cl-finite-math-only -cl-fast-relaxed-math ", storageSizeBound, storageSizePaddedBound, num_groupsBound, dims);

        if (ierr < 0) {
          printf("Error in Build Options");
          exit(-1);
        }

        err = clBuildProgram(program, 0, NULL, buildOptions, NULL, NULL);

        if (err != CL_SUCCESS) {
          std::cout << "OCL Error: compileLTwoDotBound OpenCL Build Error. Error Code: " << err << std::endl;

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

      }

      void CompileLTwoDotBoundKernels() {
        for (unsigned int i = 0; i < num_devices; ++i) {
          CompileLTwoDotBound(i, "multKernel", LTwoDotBoundKernel);

          CompileReduceBound(i, "ReduceBoundKernel", ReduceBoundKernel);
        }
      }

      void SetArgumentsLTwoDotBound() {
        cl_int ciErrNum = CL_SUCCESS;
        int counter = 0;

        for (unsigned int i = 0; i < num_devices; ++i) {
          ciErrNum |= clSetKernelArg(LTwoDotBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLevelBound[i]);
          ciErrNum |= clSetKernelArg(LTwoDotBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLevelIndexLevelintconBound[i]);

          ciErrNum |= clSetKernelArg(LTwoDotBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrIndexBound[i]);
          ciErrNum |= clSetKernelArg(LTwoDotBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLevel_intBound[i]);
          ciErrNum |= clSetKernelArg(LTwoDotBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrAlphaBound[i]);
          ciErrNum |= clSetKernelArg(LTwoDotBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrParResultBound[i]);
          ciErrNum |= clSetKernelArg(LTwoDotBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLcl_qBound[i]);
          oclCheckErr(ciErrNum, "clSetKernelArg1 Kernel Construct");

          counter = 0;
          ciErrNum |= clSetKernelArg(ReduceBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrResultBound[i]);
          ciErrNum |= clSetKernelArg(ReduceBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrParResultBound[i]);
          counter = 0;
          oclCheckErr(ciErrNum, "clSetKernelArg1 Kernel Construct");
        }

      }

      void CleanUpLTwoDotBound() {
        if (!isFirstTimeLTwoDotBound) {
          cl_int ciErrNum = CL_SUCCESS;

          for (unsigned int i = 0; i < num_devices; i++)  {
            ciErrNum |= clReleaseKernel( LTwoDotBoundKernel[i] );
            oclCheckErr(ciErrNum, "clReleaseKernel");
          }
        }
      }

      Timing PrintGFLOPSLTwoDotBound() {
        double dimsD = (double)dims;
        double StorDInnersize = (double)storageInnerSizeBound;
        double StorDBound = (double)storageSizeBound;
        double LTwoDotBoundTime = (MultTimeLTwoDotBound + ReduTimeLTwoDotBound) * (1e-9);
        LTwoDotBoundTime = LTwoDotBoundExecTime;

        if (!isFirstTimeLTwoDotBound) {
          std::cout << "Overall Runtime of AllReduceLTwoDotBound: " << LTwoDotBoundAllReduceTime << std::endl;

          std::cout << "Time for GPU-calculation of LTwoDotBound: " << LTwoDotBoundTime << std::endl;
          //        std::cout << "Number of Calls of LTwoDotBound             : " << CounterLTwoDotBound << std::endl;

        }

        double numberOfFlopBound = StorDInnersize * StorDBound * (dimsD * (28.0 + dimsD));


        if (!isFirstTimeLTwoDotBound) {
          double LTwoDotBoundGFLOPS = ((numberOfFlopBound * CounterLTwoDotBound) / LTwoDotBoundTime) * (1e-9);
          double LTwoDotBoundOPS = StorDInnersize * StorDBound * (dimsD * (20.0 + dimsD));
          double LTwoDotBoundGOPS = ((LTwoDotBoundOPS * CounterLTwoDotBound) / LTwoDotBoundTime) * (1e-9);
          return (Timing(LTwoDotBoundGFLOPS, LTwoDotBoundGOPS, LTwoDotBoundTime));
        }

        return (Timing());
      } // PrintGFLOPSLTwoDotBound
    } // namespace oclpdekernels
  }
}

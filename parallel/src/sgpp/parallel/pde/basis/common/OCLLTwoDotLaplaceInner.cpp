// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/pde/basis/common/OCLLTwoDotLaplaceInner.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    namespace oclpdekernels {

      cl_kernel LTwoDotLaplaceInnerKernel[NUMDEVS];


      double LTwoDotLaplaceInnerAllReduceTime = 0.0;
      double LTwoDotLaplaceInnerExecTime = 0.0;
      double LTwoDotLaplaceInnerReduceTime = 0.0;
      double LTwoDotLaplaceInnerExecTimeFirst = 0.0;
      double LTwoDotLaplaceInnerExecTimeLast = 0.0;
      double CounterLTwoDotLaplaceInner = 0.0;
      double CounterLTwoDotLaplaceInnerFirst = 0.0;
      double CounterLTwoDotLaplaceInnerLast = 0.0;
      double CounterLTwoDotLaplaceInnerRest = 0.0;
      double LTwoDotLaplaceInnerStartupTime = 0.0;
      double* LTwoDotLaplaceInnerAll = new double[MOD];
      double* LTwoDotLaplaceInnerProfiling  = new double[MOD];
      double* LTwoDotLaplaceInnerWaiting = new double[MOD];

      double LTwoDotLaplaceInnerProfilingAcc = 0.0;
      double LTwoDotLaplaceInnerProfilingWait = 0.0;
      void SetArgumentsLTwoDotLaplaceInner() {
        for (size_t i = 0; i < MOD; i++) {
          LTwoDotLaplaceInnerAll[i] = 0.0;
        }

        for (size_t i = 0; i < MOD; i++) {
          LTwoDotLaplaceInnerProfiling[i] = 0.0;
        }

        for (size_t i = 0; i < MOD; i++) {
          LTwoDotLaplaceInnerProfiling[i] = 0.0;
        }

        cl_int ciErrNum = CL_SUCCESS;
        int counter = 0;

        for (unsigned int i = 0; i < num_devices; ++i) {
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLevelInner[i]);
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLevelIndexLevelintconInner[i]);

          ciErrNum |= clSetKernelArg(LTwoDotLaplaceInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrIndexInner[i]);
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLevel_intInner[i]);

          ciErrNum |= clSetKernelArg(LTwoDotLaplaceInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrAlphaInner[i]);
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrParResultInner[i]);
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLcl_qInner[i]);
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLambdaInner[i]);
          oclCheckErr(ciErrNum, "clSetKernelArg1 Kernel Construct");

          counter = 0;
          ciErrNum |= clSetKernelArg(ReduceInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrResultInner[i]);
          ciErrNum |= clSetKernelArg(ReduceInnerKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrParResultInner[i]);
          counter = 0;
          oclCheckErr(ciErrNum, "clSetKernelArg1 Kernel Construct");
        }
      }

      Timing PrintGFLOPSLTwoDotLaplaceInner() {

        double StorD = (double)storageSize;
        double dimsD = (double)dims;

        double totalTime = LTwoDotLaplaceInnerExecTime;

        if (!isFirstTimeLTwoDotLaplaceInner) {
          std::cout << "Overall Runtime of AllReduceLTwoDotLaplaceInner: " << LTwoDotLaplaceInnerAllReduceTime << std::endl;
          //    std::cout << "Time for GPU-calculation of LTwoDotLaplaceInnerRedu: " << LTwoDotLaplaceInnerExecTimeFirst*(1e-9) << std::endl;
          //    std::cout << "Time for GPU-calculation of LTwoDotLaplaceInnerRedu2: " << LTwoDotLaplaceInnerExecTime*(1e-9) << std::endl;
          //    std::cout << "Time for GPU-calculation of LTwoDotLaplaceInnerProfilingAcc: " << LTwoDotLaplaceInnerProfilingAcc*(1e-9) <<  std::endl;
          //    std::cout << "Time for GPU-calculation of LTwoDotLaplaceInnerProfilingWait: " << LTwoDotLaplaceInnerProfilingWait*(1e-9) <<  std::endl;
#ifdef USE_MPI

          for (size_t i = 0; i < MOD; i++) {
            //      std::cout << "Time for GPU-calculation of LTwoDotLaplaceInner["<<i<<"]: "<< LTwoDotLaplaceInnerAll[i] << std::endl;
            //      std::cout << "   Profiling["<<i<<"]: " << LTwoDotLaplaceInnerProfiling[i]*(1e-9) ;
            //      std::cout << "   Waiting["<<i<<"]: " << LTwoDotLaplaceInnerWaiting[i] *(1e-9) << std::endl;
          }

#else

          for (size_t i = 0; i < MOD; i++) {
            // std::cout << "Time for GPU-calculation of LTwoDotLaplaceInner["<<i<<"]: "<< LTwoDotLaplaceInnerAll[i] << std::endl;
          }

#endif

          //    std::cout << "Time for GPU-calculation of LTwoDotLaplaceInnerRest:  " << LTwoDotLaplaceInnerExecTime << std::endl;
          std::cout << "Time for GPU-calculation of LTwoDotLaplaceInner: " << totalTime << std::endl;
#ifdef USE_MPI
          //        std::cout << "Number of Calls of LTwoDotLaplaceInnerFirst         : " << CounterLTwoDotLaplaceInnerFirst << std::endl;
          //        std::cout << "Number of Calls of LTwoDotLaplaceInnerRest         : " << (CounterLTwoDotLaplaceInnerRest-CounterLTwoDotLaplaceInnerLast) << std::endl;
          //        std::cout << "Number of Calls of LTwoDotLaplaceInnerLast         : " << CounterLTwoDotLaplaceInnerLast << std::endl;

          //        std::cout << "Number of Calls of LTwoDotLaplaceInner             : " << CounterLTwoDotLaplaceInner << std::endl;

          //    std::cout << "Time for Single GPU-calculation of LTwoDotLaplaceInnerFirst: " << LTwoDotLaplaceInnerExecTimeFirst/CounterLTwoDotLaplaceInnerFirst << std::endl;
          //    std::cout << "Time for Single GPU-calculation of LTwoDotLaplaceInnerRest : " << LTwoDotLaplaceInnerExecTime/(CounterLTwoDotLaplaceInnerRest-CounterLTwoDotLaplaceInnerLast) << std::endl;
          //    std::cout << "Time for Single GPU-calculation of LTwoDotLaplaceInnerLast : " << LTwoDotLaplaceInnerExecTimeLast/CounterLTwoDotLaplaceInnerLast << std::endl;
#endif
        }

        // LTwoDotLaplace Inner
        double numberOfFlop = StorD * StorD * (dimsD * (27.0 + dimsD) + 2);

        if (!isFirstTimeLTwoDotLaplaceInner) {
          double LTwoDotLaplaceInnerGFLOPS = ((numberOfFlop * CounterLTwoDotLaplaceInner) / totalTime) * (1e-9);
          double LTwoDotLaplaceInnerOPS = StorD * StorD * (dimsD * (16.0));
          double LTwoDotLaplaceInnerGOPS = ((LTwoDotLaplaceInnerOPS * CounterLTwoDotLaplaceInner) / totalTime) * (1e-9);

          //    std::cout << "GFLOPS for LTwoDotLaplaceInner : " << LTwoDotLaplaceInnerGFLOPS << std::endl;

          return (Timing(LTwoDotLaplaceInnerGFLOPS, LTwoDotLaplaceInnerGOPS, totalTime));
        }

        return (Timing());
      }


      std::string LTwoDotLaplaceInnerHeader() {
        std::stringstream stream_program_src;
        stream_program_src << "__kernel void multKernel(__global  REAL* ptrLevel,			"  << std::endl;
        stream_program_src << "			 	__constant  REAL* ptrLevelIndexLevelintcon,	"  << std::endl;
        stream_program_src << "			 	__global  REAL* ptrIndex, 			"  << std::endl;
        stream_program_src << "			 	__global  REAL* ptrLevel_int, 			"  << std::endl;
        stream_program_src << "			 	__global  REAL* ptrAlpha, 			"  << std::endl;
        stream_program_src << "			 	__global  REAL* ptrParResult, 			"  << std::endl;
        stream_program_src << "			 	__constant  REAL* ptrLcl_q,			"  << std::endl;
        stream_program_src << "			 	__constant  REAL* ptrLambda,			"  << std::endl;
        stream_program_src << "			 	ulong overallMultOffset,   			"  << std::endl;
        stream_program_src << "                        	ulong ConstantMemoryOffset,   			"  << std::endl;
        stream_program_src << "                        	REAL TimestepCoeff)		                "  << std::endl;
        stream_program_src << "{									"  << std::endl;
        return stream_program_src.str();
      }



      void CompileLTwoDotLaplaceInner(int id, std::string kernel_src, cl_kernel* kernel) {
        cl_int err = CL_SUCCESS;

        std::string source2 = LTwoDotLaplaceInnerHeader();
        std::string l2dotfunction = InnerLTwoDotFunction();
        std::string gradientfunction = InnerGradientFunction();
        std::stringstream stream_program_src;


#ifdef USEOCL_NVIDIA
        stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
#endif
        stream_program_src << "#define REAL double" << std::endl;
        stream_program_src << "#define LSIZE " << LSIZE << std::endl;
        stream_program_src << l2dotfunction << std::endl;
        stream_program_src << gradientfunction << std::endl;
        stream_program_src << source2 << std::endl;
        stream_program_src << "__local REAL alphaTemp[" << LSIZE << "];" << std::endl;
        cl_uint sizec;
        clGetDeviceInfo(device_ids[0], CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, sizeof(cl_uint), &sizec, 0);

        size_t keplerUseRegister = (sizec >= 3);

        size_t fermiUseRegister = (dims <= 5);

#ifdef USEOCL_INTEL
        keplerUseRegister = 1;
        fermiUseRegister = 0;
#endif
        size_t dcon = (fermiUseRegister || keplerUseRegister) ? 0 : dims;

        stream_program_src << "__local REAL l2dotTemp[" << dcon* LSIZE << "];" << std::endl;
        stream_program_src << "__local REAL gradTemp[" << dcon* LSIZE << "];" << std::endl;

        stream_program_src << "ptrLevel += get_local_id(0);" << std::endl;
        stream_program_src << "ptrIndex += get_local_id(0);" << std::endl;
        stream_program_src << "ptrLevel_int += get_local_id(0);" << std::endl;

        stream_program_src << "alphaTemp[get_local_id(0)]   = ptrAlpha[get_global_id(0) + overallMultOffset + ConstantMemoryOffset*LSIZE];" << std::endl;

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

          if (d_inner < dcon) {
            stream_program_src << "l2dotTemp[" << d_inner* LSIZE << " + get_local_id(0)] = (l2dot(i_level" << registerIdStr << ",\
             i_index" << registerIdStr << ",\
             i_level_int" << registerIdStr << ", \
            j_levelreg,j_indexreg,j_level_intreg," << "ptrLcl_q[" << (d_inner + 1) * 2 - 2 << "]));" << std::endl;

            stream_program_src << "gradTemp[" << d_inner* LSIZE << " + get_local_id(0)] = (gradient(i_level" << registerIdStr << ",\
               i_index" << registerIdStr << ", \
               j_levelreg,j_indexreg," << "ptrLcl_q[" << (d_inner + 1) * 2 - 1 << "]));" << std::endl;
          } else {
            stream_program_src << "REAL l2dotreg_" << d_inner << " = \
	                 (l2dot(i_level" << registerIdStr << ",\
	                        i_index" << registerIdStr << ",\
	                        i_level_int" << registerIdStr << ",\
	                        j_levelreg,j_indexreg,j_level_intreg," << "ptrLcl_q[" << (d_inner + 1) * 2 - 2 << "]));" << std::endl;

            stream_program_src << "REAL gradreg_" << d_inner << " = \
	               (gradient(i_level" << registerIdStr << ",\
	                         i_index" << registerIdStr << ", \
	                         j_levelreg,j_indexreg," << "ptrLcl_q[" << (d_inner + 1) * 2 - 1 << "]));" << std::endl;
          }

        }

        stream_program_src << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
        typeREAL = "REAL ";
        stream_program_src << typeREAL << "alphaReg = alphaTemp[" << "k" << "];" << std::endl;

        for (unsigned d_outer = 0; d_outer < dims; d_outer++) {
          if (d_outer == 0) {
            typeREAL = "REAL ";
          } else {
            typeREAL = "";
          }

          stream_program_src << typeREAL << "element = alphaReg;" << std::endl;

          for (unsigned d_inner = 0; d_inner < dims; d_inner++) {
            if (d_outer == d_inner) {

            } else {
              stream_program_src << "element *= ";

              if (d_inner < dcon) {
                stream_program_src << "(l2dotTemp[" << d_inner* LSIZE << " + get_local_id(0)]);";
              } else {
                stream_program_src << "(l2dotreg_" << d_inner << ");";
              }
            }

            stream_program_src << std::endl;
          }

          if (d_outer < dcon) {

            stream_program_src << "res += TimestepCoeff*ptrLambda[" << d_outer << "] * element* gradTemp[" << d_outer* LSIZE << " + get_local_id(0)];" << std::endl;
          } else {

            stream_program_src << "res += TimestepCoeff*ptrLambda[" << d_outer << "] * element* gradreg_" << d_outer << ";" << std::endl;
          }

        }

        if ((dims - 1) < dcon) {
          stream_program_src << "res += element * l2dotTemp[" << (dims - 1)*LSIZE << " + get_local_id(0)]; " << std::endl;
        } else {
          stream_program_src << "res += element * (l2dotreg_" << (dims - 1) << "); " << std::endl;
        }

        stream_program_src << "}" << std::endl;
        stream_program_src << "ptrParResult[((get_group_id(0)+ConstantMemoryOffset)*get_global_size(1) + get_global_id(1))*" << LSIZE << " + get_local_id(0)] = res; " << std::endl;
        stream_program_src << "}" << std::endl;

        std::string program_src = stream_program_src.str();
        const char* source3 = program_src.c_str();
#if PRINTOCL
        std::cout <<  source3 << std::endl;
#endif
        cl_program program = clCreateProgramWithSource(context, 1, (const char**)&source3, NULL, &err);
        oclCheckErr(err, "clCreateProgramWithSource");
        char buildOptions[256];
        int ierr = snprintf(buildOptions, sizeof(buildOptions), "-DSTORAGE=%lu -DSTORAGEPAD=%lu -DNUMGROUPS=%lu -DDIMS=%lu -cl-finite-math-only -cl-fast-relaxed-math ", storageSize, storageSizePadded, num_groups, dims);

        if (ierr < 0) {
          printf("Error in Build Options");
          exit(-1);
        }

        err = clBuildProgram(program, 0, NULL, buildOptions, NULL, NULL);

        if (err != CL_SUCCESS) {
          std::cout << "OCL Error: compileLTwoDotLaplaceinner OpenCL Build Error. Error Code: " << err << std::endl;

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


      void CompileLTwoDotLaplaceInnerKernels() {
        for (unsigned int i = 0; i < num_devices; ++i) {
          CompileLTwoDotLaplaceInner(i, "multKernel", LTwoDotLaplaceInnerKernel);
          CompileReduceInner(i, "ReduceInnerKernel", ReduceInnerKernel);
        }
      }

      void CleanUpLTwoDotLaplaceInner() {
        if (!isFirstTimeLTwoDotLaplaceInner) {
          cl_int ciErrNum = CL_SUCCESS;

          for (unsigned int i = 0; i < num_devices; i++)  {
            ciErrNum |= clReleaseKernel( LTwoDotLaplaceInnerKernel[i] );
            oclCheckErr(ciErrNum, "clReleaseMemObject & clReleaseKernel");
          }

          if (!isFirstTimeLaplaceInner) {
            for (unsigned int i = 0; i < num_devices; i++)  {
              ciErrNum |= clReleaseMemObject(d_ptrLambdaInner[i]);
            }
          }

        }
      } // CleanUpLTwoDotLaplaceInner


    }
  }
}


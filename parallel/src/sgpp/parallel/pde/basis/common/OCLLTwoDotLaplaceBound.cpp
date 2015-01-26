// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include "OCLLTwoDotLaplaceBound.hpp"

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    namespace oclpdekernels {

      cl_kernel LTwoDotLaplaceBoundKernel[NUMDEVS];


      double LTwoDotLaplaceBoundAllReduceTime = 0.0;
      double LTwoDotLaplaceBoundExecTime = 0.0;
      double CounterLTwoDotLaplaceBound = 0.0;
      double LTwoDotLaplaceBoundStartupTime = 0.0;
      double* LTwoDotLaplaceBoundAll = new double[MOD];


      void SetArgumentsLTwoDotLaplaceBound() {

        for (size_t i = 0; i < MOD; i++) {
          LTwoDotLaplaceBoundAll[i] = 0.0;
        }


        cl_int ciErrNum = CL_SUCCESS;
        int counter = 0;

        for (unsigned int i = 0; i < num_devices; ++i) {
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLevelBound[i]);
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLevelIndexLevelintconBound[i]);

          ciErrNum |= clSetKernelArg(LTwoDotLaplaceBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrIndexBound[i]);
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLevel_intBound[i]);
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrParResultBound[i]);
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrAlphaBound[i]);
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLcl_qBound[i]);
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrLambdaBound[i]);
          oclCheckErr(ciErrNum, "clSetKernelArg1 LTwoDotLaplaceBoundKernel Construct");

          counter = 0;
          ciErrNum |= clSetKernelArg(ReduceBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrResultBound[i]);
          ciErrNum |= clSetKernelArg(ReduceBoundKernel[i], counter++, sizeof(cl_mem), (void*) &d_ptrParResultBound[i]);
          counter = 0;
          oclCheckErr(ciErrNum, "clSetKernelArg1 ReduceKernel Construct");
        }
      }

      Timing PrintGFLOPSLTwoDotLaplaceBound() {

        double dimsD = (double)dims;
        double StorDInnersize = (double)storageInnerSizeBound;
        double StorDBound = (double)storageSizeBound;

        if (!isFirstTimeLTwoDotLaplaceBound) {
          std::cout << "Overall Runtime of AllReduceLTwoDotLaplaceBound: " << LTwoDotLaplaceBoundAllReduceTime << std::endl;
          std::cout << "Time for GPU-calculation of LTwoDotLaplaceBound: " << LTwoDotLaplaceBoundExecTime << std::endl;
          //        std::cout << "Number of Calls of LTwoDotLaplaceBound             : " << CounterLTwoDotLaplaceBound << std::endl;
#ifdef USE_MPI

          for (size_t i = 0; i < MOD; i++) {
            //      std::cout << "Time for GPU-calculation of LTwoDotLaplaceBound["<<i<<"]: "<< LTwoDotLaplaceBoundAll[i] << std::endl;
          }

#endif
        }

        // LTwoDotLaplace Bound
        double numberOfFlop = StorDInnersize * StorDBound * (dimsD * (29.0 + dimsD) + 2);

        if (!isFirstTimeLTwoDotLaplaceBound) {
          double LTwoDotLaplaceBoundGFLOPS = ((numberOfFlop * CounterLTwoDotLaplaceBound) / LTwoDotLaplaceBoundExecTime) * (1e-9);
          double LTwoDotLaplaceBoundOPS = StorDInnersize * StorDBound * (dimsD * (26.0));
          double LTwoDotLaplaceBoundGOPS = ((LTwoDotLaplaceBoundOPS * CounterLTwoDotLaplaceBound) / LTwoDotLaplaceBoundExecTime) * (1e-9);
          //    std::cout << "GFLOPS for LTwoDotLaplaceBound : " << LTwoDotLaplaceBoundGFLOPS << std::endl;
          return (Timing(LTwoDotLaplaceBoundGFLOPS, LTwoDotLaplaceBoundGOPS, LTwoDotLaplaceBoundExecTime));
        }

        return (Timing());
      }


      std::string LTwoDotLaplaceBoundHeader() {
        std::stringstream stream_program_src;
        stream_program_src << "__kernel void multKernel(	__global  REAL* ptrLevel,       "  << std::endl;
        stream_program_src << "			 		__constant  REAL* ptrLevelIndexLevelintcon,"  << std::endl;
        stream_program_src << "			 		__global  REAL* ptrIndex, 	"  << std::endl;
        stream_program_src << "					__global  REAL* ptrLevel_int, 	"  << std::endl;
        stream_program_src << "				 	__global  REAL* ptrParResult, 	"  << std::endl;
        stream_program_src << "			 		__global  REAL* ptrAlpha, 	"  << std::endl;
        stream_program_src << "					__constant  REAL* ptrLcl_q,	"  << std::endl;
        stream_program_src << "				 	__constant  REAL* ptrLambda, 	"  << std::endl;
        stream_program_src << "	    				ulong overallMultOffset,	"  << std::endl;
        stream_program_src << "				 	ulong ConstantMemoryOffset,	"  << std::endl;
        stream_program_src << "                        	REAL TimestepCoeff)		                "  << std::endl;
        stream_program_src << "{									"  << std::endl;
        return stream_program_src.str();
      }




      void CompileLTwoDotLaplaceBound(int id, std::string kernel_src, cl_kernel* kernel) {


        cl_int err = CL_SUCCESS;
        std::string source2 = LTwoDotLaplaceBoundHeader();
        std::string l2dotfunction = BoundLTwoDotFunction();
        std::string gradientfunction = BoundGradientFunction();

        std::stringstream stream_program_src;
#ifdef USEOCL_NVIDIA
        stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
#endif
        stream_program_src << "#define REAL double" << std::endl;
        stream_program_src << "#define LSIZE " << LSIZE << std::endl;
        stream_program_src << l2dotfunction << std::endl;
        stream_program_src << gradientfunction << std::endl;
        stream_program_src << source2 << std::endl;
        stream_program_src << "	 __local REAL alphaTemp[" << LSIZE << "];	"  << std::endl;

        //CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV
        cl_uint sizec;
        clGetDeviceInfo(device_ids[0], CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, sizeof(cl_uint), &sizec, 0);

        size_t keplerUseRegister = (sizec >= 3);
        size_t fermiUseRegister = (dims <= 5);

#ifdef USEOCL_INTEL
        keplerUseRegister = 0;
        fermiUseRegister = 0;
#endif

        size_t dcon = (fermiUseRegister || keplerUseRegister) ? 0 : dims;

        stream_program_src << "__local REAL l2dotTemp[" << dcon* LSIZE << "];" << std::endl;
        stream_program_src << "__local REAL gradTemp[" << dcon* LSIZE << "];" << std::endl;


        stream_program_src << "ptrLevel += get_local_id(0);" << std::endl;
        stream_program_src << "ptrIndex += get_local_id(0);" << std::endl;
        stream_program_src << "ptrLevel_int += get_local_id(0);" << std::endl;

        stream_program_src << "alphaTemp[get_local_id(0)]   = ptrAlpha[get_global_id(0) +overallMultOffset +  ConstantMemoryOffset*LSIZE];" << std::endl;

        stream_program_src << "REAL res = 0.0;" << std::endl;

        if (keplerUseRegister) {
          const char* typeREAL2 = "REAL ";

          for (size_t d_inner = 0; d_inner < dims; d_inner++) {
            stream_program_src << typeREAL2 << "i_level_" << d_inner << " =         ptrLevel[" << d_inner* storageInnerSizePaddedBound << "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
            stream_program_src << typeREAL2 << "i_index_" << d_inner << " =         ptrIndex[" << d_inner* storageInnerSizePaddedBound << "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
            stream_program_src << typeREAL2 << "i_level_int_" << d_inner << " =     ptrLevel_int[" << d_inner* storageInnerSizePaddedBound << "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
          }
        }


        const char* typeREAL = "";

        stream_program_src << "for (unsigned k = 0; k < " << LSIZE << "; k++) {" << std::endl;

        for (size_t d_inner = 0; d_inner < dims; d_inner++) {
          if (d_inner == 0) {
            typeREAL = "REAL ";
          } else {
            typeREAL = "";
          }

          if (!keplerUseRegister) {
            const char* typeREAL2 = "REAL ";
            stream_program_src << typeREAL2 << "i_level_" << d_inner << " =         ptrLevel[" << d_inner* storageInnerSizePaddedBound << "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
            stream_program_src << typeREAL2 << "i_index_" << d_inner << " =         ptrIndex[" << d_inner* storageInnerSizePaddedBound << "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
            stream_program_src << typeREAL2 << "i_level_int_" << d_inner << " =     ptrLevel_int[" << d_inner* storageInnerSizePaddedBound << "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
          }

          stream_program_src << typeREAL << "j_levelreg = ptrLevelIndexLevelintcon[(get_group_id(0)*LSIZE+k)*" << dims * 3 << " + " << d_inner * 3 << "]; " << std::endl;
          stream_program_src << typeREAL << "j_indexreg = ptrLevelIndexLevelintcon[(get_group_id(0)*LSIZE+k)*" << dims * 3 << " + " << d_inner * 3 + 1 << "]; " << std::endl;
          stream_program_src << typeREAL << "j_level_intreg = ptrLevelIndexLevelintcon[(get_group_id(0)*LSIZE+k)*" << dims * 3 << " + " << d_inner * 3 + 2 << "]; " << std::endl;


          if (d_inner < dcon) {
            stream_program_src << "l2dotTemp[" << d_inner* LSIZE << " + get_local_id(0)] = (l2dot(i_level_" << d_inner << ",\
             i_index_" << d_inner << ",\
             i_level_int_" << d_inner << ",\
             j_levelreg,j_indexreg,j_level_intreg," << "ptrLcl_q[" << (d_inner + 1) * 2 - 2 << "]));" << std::endl;

            stream_program_src << "gradTemp[" << d_inner* LSIZE << " + get_local_id(0)] = (gradient(i_level_" << d_inner << ",\
               i_index_" << d_inner << ",\
               j_levelreg,j_indexreg," << "ptrLcl_q[" << (d_inner + 1) * 2 - 1 << "]));" << std::endl;
          } else {
            stream_program_src << "REAL l2dotreg_" << d_inner << " = \
            (l2dot(i_level_" << d_inner << ",\
                   i_index_" << d_inner << ",\
                   i_level_int_" << d_inner << ",\
                   j_levelreg,j_indexreg,j_level_intreg," << "ptrLcl_q[" << (d_inner + 1) * 2 - 2 << "]));" << std::endl;

            stream_program_src << "REAL gradreg_" << d_inner << " = \
            (gradient(i_level_" << d_inner << ",\
                      i_index_" << d_inner << ",\
                      j_levelreg,j_indexreg," << "ptrLcl_q[" << (d_inner + 1) * 2 - 1 << "]));" << std::endl;
          }

        }

        stream_program_src << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
        typeREAL = "REAL ";
        stream_program_src << typeREAL << "alphaTempReg = alphaTemp[" << "k" << "];" << std::endl;

        for (unsigned d_outer = 0; d_outer < dims; d_outer++) {
          if (d_outer == 0) {
            typeREAL = "REAL ";
          } else {
            typeREAL = "";
          }

          stream_program_src << typeREAL << "element = alphaTempReg;" << std::endl;

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
          std::cout << "OCL Error: compileLTwoDotLaplaceBoundKernel OpenCL Build Error. Error Code: " << err << std::endl;

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

      } // compileLTwoDotLaplaceBoundKernel


      void CompileLTwoDotLaplaceBoundKernels() {
        for (unsigned int i = 0; i < num_devices; ++i) {
          CompileLTwoDotLaplaceBound(i, "multKernel", LTwoDotLaplaceBoundKernel);
          CompileReduceBound(i, "ReduceBoundKernel", ReduceBoundKernel);
        }
      }

      void CleanUpLTwoDotLaplaceBound() {
        if (!isFirstTimeLTwoDotLaplaceBound) {
          cl_int ciErrNum = CL_SUCCESS;

          for (unsigned int i = 0; i < num_devices; i++)  {
            ciErrNum |= clReleaseKernel( LTwoDotLaplaceBoundKernel[i] );
            oclCheckErr(ciErrNum, "clReleaseMemObject & clReleaseKernel");
          }

          if (!isFirstTimeLaplaceBound) {
            for (unsigned int i = 0; i < num_devices; i++)  {
              ciErrNum |= clReleaseMemObject(d_ptrLambdaBound[i]);
            }
          }

        }
      }


    }
  }
}


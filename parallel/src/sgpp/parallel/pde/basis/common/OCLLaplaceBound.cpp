// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OCLPDEBound.hpp"

using namespace SGPP::base;
#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace parallel {
namespace oclpdekernels {

cl_kernel LaplaceBoundKernel[NUMDEVS];
cl_mem d_ptrLambdaBound[NUMDEVS];
double MultTimeLaplaceBound = 0.0;
double ReduTimeLaplaceBound = 0.0;
double CounterLaplaceBound = 0.0;
double LaplaceBoundStartupTime = 0.0;
double LaplaceBoundExecTime = 0.0;
double LaplaceBoundAllReduceTime = 0.0;
//SGppStopwatch* myStopwatch = new SGppStopwatch();


void SetLambdaBufferLaplaceBound(REAL* ptrLambda,
                                 size_t localdim) {
  size_t lambda_size = localdim;
  cl_int ciErrNum = CL_SUCCESS;

  for (unsigned int i = 0; i < num_devices; ++i) {

    d_ptrLambdaBound[i] = clCreateBuffer(context,
                                         CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                         lambda_size * sizeof(REAL), ptrLambda, &ciErrNum);
    oclCheckErr(ciErrNum, "clCreateBuffer");

  }

}

std::string BoundGradientFunction() {
  std::stringstream stream_program_src;
  stream_program_src << "REAL gradient(	REAL i_level, 	"  << std::endl;
  stream_program_src << "			REAL i_index, 	"  << std::endl;
  stream_program_src << "			REAL j_level, 	"  << std::endl;
  stream_program_src << "			REAL j_index,	"  << std::endl;
  stream_program_src << "			REAL lcl_q_inv)	"  << std::endl;
  stream_program_src << "{				"  << std::endl;
  stream_program_src << "  REAL grad;			"  << std::endl;
  stream_program_src <<
                     "  ulong doGrad = (ulong)((i_level == j_level) && (i_index == j_index) && (i_level != 1.0));"
                     << std::endl;
  stream_program_src <<
                     "  grad = select(0.0, i_level * 2.0 * lcl_q_inv, doGrad);		"  << std::endl;
  stream_program_src << "  return grad;							"  << std::endl;
  stream_program_src << "}								"  << std::endl;
  return stream_program_src.str();
}

std::string LaplaceBoundHeader() {
  std::stringstream stream_program_src;
  stream_program_src <<
                     "__kernel void multKernelBound(	__global  REAL* ptrLevel,       "  << std::endl;
  stream_program_src << "			 		__constant  REAL* ptrLevelIndexLevelintcon,"  <<
                     std::endl;
  stream_program_src << "			 		__global  REAL* ptrIndex, 	"  << std::endl;
  stream_program_src << "					__global  REAL* ptrLevel_int, 	"  << std::endl;
  stream_program_src << "				 	__global  REAL* ptrParResult, 	"  << std::endl;
  stream_program_src << "			 		__global  REAL* ptrAlpha, 	"  << std::endl;
  stream_program_src << "					__constant  REAL* ptrLcl_q,	"  << std::endl;
  stream_program_src << "				 	__constant  REAL* ptrLambda, 	"  << std::endl;
  stream_program_src << "	    				ulong overallMultOffset,	"  << std::endl;
  stream_program_src << "				 	ulong ConstantMemoryOffset)	"  << std::endl;
  stream_program_src << "{								"  << std::endl;
  return stream_program_src.str();
}


void CompileLaplaceBound(int id, std::string kernel_src, cl_kernel* kernel) {
  cl_int err = CL_SUCCESS;
  std::string source2 = LaplaceBoundHeader();
  std::string l2dotfunction = BoundLTwoDotFunction();
  std::string gradientfunction = BoundGradientFunction();

  std::stringstream stream_program_src;
#ifdef USEOCL_NVIDIA
  stream_program_src << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" <<
                     std::endl << std::endl;
#endif
  stream_program_src << "#define REAL double" << std::endl;
  stream_program_src << "#define LSIZE " << LSIZE << std::endl;
  stream_program_src << l2dotfunction << std::endl;
  stream_program_src << gradientfunction << std::endl;
  stream_program_src << source2 << std::endl;
  stream_program_src << "	 __local REAL alphaTemp[LSIZE];	"  << std::endl;


  cl_uint sizec;
  clGetDeviceInfo(device_ids[0], CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV,
                  sizeof(cl_uint), &sizec, 0);

  size_t keplerUseRegister = (sizec >= 3);
  size_t fermiUseRegister = (dims <= 5);

#ifdef USEOCL_INTEL
  keplerUseRegister = 0;
  fermiUseRegister = 0;
#endif

  size_t dcon = (fermiUseRegister || keplerUseRegister) ? 0 : dims;



  stream_program_src << "__local REAL l2dotTemp[" << dcon* LSIZE << "];" <<
                     std::endl;
  stream_program_src << "__local REAL gradTemp[" << dcon* LSIZE << "];" <<
                     std::endl;


  stream_program_src << "ptrLevel += get_local_id(0);" << std::endl;
  stream_program_src << "ptrIndex += get_local_id(0);" << std::endl;
  stream_program_src << "ptrLevel_int += get_local_id(0);" << std::endl;

  stream_program_src <<
                     "alphaTemp[get_local_id(0)]   = ptrAlpha[get_global_id(0) +overallMultOffset +  ConstantMemoryOffset*LSIZE];"
                     << std::endl;

  stream_program_src << "REAL res = 0.0;" << std::endl;

  if (keplerUseRegister) {
    const char* typeREAL2 = "REAL ";

    for (size_t d_inner = 0; d_inner < dims; d_inner++) {
      stream_program_src << typeREAL2 << "i_level_" << d_inner <<
                         " =         ptrLevel[" << d_inner* storageInnerSizePaddedBound <<
                         "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
      stream_program_src << typeREAL2 << "i_index_" << d_inner <<
                         " =         ptrIndex[" << d_inner* storageInnerSizePaddedBound <<
                         "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
      stream_program_src << typeREAL2 << "i_level_int_" << d_inner <<
                         " =     ptrLevel_int[" << d_inner* storageInnerSizePaddedBound <<
                         "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
    }
  }

  const char* typeREAL = "";

  stream_program_src << "for (unsigned k = 0; k < " << LSIZE << "; k++) {" <<
                     std::endl;

  for (size_t d_inner = 0; d_inner < dims; d_inner++) {
    if (d_inner == 0) {
      typeREAL = "REAL ";
    } else {
      typeREAL = "";
    }

    if (!keplerUseRegister) {
      const char* typeREAL2 = "REAL ";
      stream_program_src << typeREAL2 << "i_level_" << d_inner <<
                         " =         ptrLevel[" << d_inner* storageInnerSizePaddedBound <<
                         "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
      stream_program_src << typeREAL2 << "i_index_" << d_inner <<
                         " =         ptrIndex[" << d_inner* storageInnerSizePaddedBound <<
                         "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
      stream_program_src << typeREAL2 << "i_level_int_" << d_inner <<
                         " =     ptrLevel_int[" << d_inner* storageInnerSizePaddedBound <<
                         "+ (get_global_id(1))*" << LSIZE << "];" << std::endl;
    }


    stream_program_src << typeREAL <<
                       "j_levelreg = ptrLevelIndexLevelintcon[(get_group_id(0)*LSIZE+k)*" << dims * 3
                       << " + " << d_inner * 3 << "]; " << std::endl;
    stream_program_src << typeREAL <<
                       "j_indexreg = ptrLevelIndexLevelintcon[(get_group_id(0)*LSIZE+k)*" << dims * 3
                       << " + " << d_inner * 3 + 1 << "]; " << std::endl;
    stream_program_src << typeREAL <<
                       "j_level_intreg = ptrLevelIndexLevelintcon[(get_group_id(0)*LSIZE+k)*" << dims *
                       3 << " + " << d_inner * 3 + 2 << "]; " << std::endl;


    if (d_inner < dcon) {
      stream_program_src << "l2dotTemp[" << d_inner* LSIZE <<
                         " + get_local_id(0)] = (l2dot(i_level_" << d_inner << ",\
             i_index_" << d_inner << ",\
             i_level_int_" << d_inner << ",\
             j_levelreg,j_indexreg,j_level_intreg," << "ptrLcl_q[" <<
                         (d_inner + 1) * 2 - 2 << "]));" << std::endl;

      stream_program_src << "gradTemp[" << d_inner* LSIZE <<
                         " + get_local_id(0)] = (gradient(i_level_" << d_inner << ",\
               i_index_" << d_inner << ",\
               j_levelreg,j_indexreg," << "ptrLcl_q[" <<
                         (d_inner + 1) * 2 - 1 << "]));" << std::endl;
    } else {
      stream_program_src << "REAL l2dotreg_" << d_inner << " = \
            (l2dot(i_level_" << d_inner << ",\
                   i_index_" << d_inner << ",\
                   i_level_int_" << d_inner << ",\
                   j_levelreg,j_indexreg,j_level_intreg," << "ptrLcl_q[" <<
                         (d_inner + 1) * 2 - 2 << "]));" << std::endl;

      stream_program_src << "REAL gradreg_" << d_inner << " = \
            (gradient(i_level_" << d_inner << ",\
                      i_index_" << d_inner << ",\
                      j_levelreg,j_indexreg," << "ptrLcl_q[" <<
                         (d_inner + 1) * 2 - 1 << "]));" << std::endl;
    }

  }

  stream_program_src << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
  typeREAL = "REAL ";
  stream_program_src << typeREAL << "alphaTempReg = alphaTemp[" << "k" << "];" <<
                     std::endl;

  for (unsigned d_outer = 0; d_outer < dims; d_outer++) {
    if (d_outer == 0) {
      typeREAL = "REAL ";
    } else {
      typeREAL = "";
    }

    stream_program_src << typeREAL << "element = alphaTempReg;" << std::endl;

    for (unsigned d_inner = 0; d_inner < dims; d_inner++) {
      stream_program_src << "element *= ";

      if (d_outer == d_inner) {
        if (d_inner < dcon) {
          stream_program_src << "(gradTemp[" << d_inner* LSIZE << " + get_local_id(0)]);";
        } else {
          stream_program_src << "(gradreg_" << d_inner << ");";
        }
      } else {
        if (d_inner < dcon) {
          stream_program_src << "(l2dotTemp[" << d_inner* LSIZE <<
                             " + get_local_id(0)]);";
        } else {
          stream_program_src << "(l2dotreg_" << d_inner << ");";
        }
      }

      stream_program_src << std::endl;
    }

    stream_program_src << "res += ptrLambda[" << d_outer << "] * element;" <<
                       std::endl;
  }

  stream_program_src << "}" << std::endl;
  stream_program_src <<
                     "ptrParResult[((get_group_id(0) + ConstantMemoryOffset)*get_global_size(1) + get_global_id(1))*"
                     << LSIZE << " + get_local_id(0)] = res; " << std::endl;
  stream_program_src << "}" << std::endl;

  std::string program_src = stream_program_src.str();
  source2 = program_src.c_str();
#if PRINTOCL
  std::cout <<  source2 << std::endl;
#endif
  cl_program program = clCreateProgramWithSource(context, 1,
                       (const char**)&source2, NULL, &err);
  oclCheckErr(err, "clCreateProgramWithSource");
  char buildOptions[256];
  int ierr = snprintf(buildOptions, sizeof(buildOptions),
                      "-DSTORAGE=%lu -DSTORAGEPAD=%lu -DNUMGROUPS=%lu -DDIMS=%lu -cl-finite-math-only -cl-fast-relaxed-math ",
                      storageSizeBound, storageSizePaddedBound, num_groupsBound, dims);

  if (ierr < 0) {
    printf("Error in Build Options");
    exit(-1);
  }

  err = clBuildProgram(program, 0, NULL, buildOptions, NULL, NULL);

  if (err != CL_SUCCESS) {
    std::cout <<
              "OCL Error: compileMultKernel2DBound2Indexes OpenCL Build Error. Error Code: "
              << err << std::endl;

    size_t len;
    char buffer[10000];

    // get the build log
    clGetProgramBuildInfo(program, device_ids[0], CL_PROGRAM_BUILD_LOG,
                          sizeof(buffer), buffer, &len);

    std::cout << "--- Build Log ---" << std::endl << buffer << std::endl;
  }

  oclCheckErr(err, "clBuildProgram");


  kernel[id] = clCreateKernel(program, kernel_src.c_str(), &err);
  oclCheckErr(err, "clCreateKernel");

  err |= clReleaseProgram(program);
  oclCheckErr(err, "clReleaseProgram");

} // compileMultKernel2DBound2Indexes

void CompileLaplaceBoundKernels() {
  for (unsigned int i = 0; i < num_devices; ++i) {
    CompileLaplaceBound(i, "multKernelBound", LaplaceBoundKernel);

    CompileReduceBound(i, "ReduceBoundKernel", ReduceBoundKernel);
  }
}

void SetArgumentsLaplaceBound() {
  cl_int ciErrNum = CL_SUCCESS;
  int counter = 0;

  for (unsigned int i = 0; i < num_devices; ++i) {
    ciErrNum |= clSetKernelArg(LaplaceBoundKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrLevelBound[i]);
    ciErrNum |= clSetKernelArg(LaplaceBoundKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrLevelIndexLevelintconBound[i]);

    ciErrNum |= clSetKernelArg(LaplaceBoundKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrIndexBound[i]);
    ciErrNum |= clSetKernelArg(LaplaceBoundKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrLevel_intBound[i]);
    ciErrNum |= clSetKernelArg(LaplaceBoundKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrParResultBound[i]);
    ciErrNum |= clSetKernelArg(LaplaceBoundKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrAlphaBound[i]);
    ciErrNum |= clSetKernelArg(LaplaceBoundKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrLcl_qBound[i]);
    ciErrNum |= clSetKernelArg(LaplaceBoundKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrLambdaBound[i]);
    oclCheckErr(ciErrNum, "clSetKernelArg1 Kernel Construct");

    counter = 0;
    ciErrNum |= clSetKernelArg(ReduceBoundKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrResultBound[i]);
    ciErrNum |= clSetKernelArg(ReduceBoundKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrParResultBound[i]);
    counter = 0;
  }
}


void CleanUpLaplaceBound() {
  if (!isFirstTimeLaplaceBound) {
    cl_int ciErrNum = CL_SUCCESS;

    for (unsigned int i = 0; i < num_devices; i++)  {
      ciErrNum |= clReleaseMemObject(d_ptrLambdaBound[i]);
      ciErrNum |= clReleaseKernel( LaplaceBoundKernel[i] );
      oclCheckErr(ciErrNum, "clReleaseMemObject & clReleaseKernel");
    }
  }
}
Timing PrintGFLOPSLaplaceBound() {
  double dimsD = (double)dims;
  double StorDInnersize = (double)storageInnerSizeBound;
  double StorDBound = (double)storageSizeBound;
  double LaplaceBoundTime = (MultTimeLaplaceBound + ReduTimeLaplaceBound) *
                            (1e-9);
  LaplaceBoundTime = LaplaceBoundExecTime;

  if (!isFirstTimeLaplaceBound) {
    std::cout << "Overall Runtime of AllReduceLaplaceBound: " <<
              LaplaceBoundAllReduceTime << std::endl;
    std::cout << "Time for GPU-calculation of LaplaceBound: " << LaplaceBoundTime <<
              std::endl;
  }

  double numberOfFlopBound = StorDInnersize * StorDBound * (dimsD *
                             (28.0 + dimsD));


  if (!isFirstTimeLaplaceBound) {
    double LaplaceBoundGFLOPS = ((numberOfFlopBound * CounterLaplaceBound) /
                                 LaplaceBoundTime) * (1e-9);
    double LaplaceBoundOPS = StorDInnersize * StorDBound * (dimsD * (26.0));
    double LaplaceBoundGOPS = ((LaplaceBoundOPS * CounterLaplaceBound) /
                               LaplaceBoundTime) * (1e-9);
    return (Timing(LaplaceBoundGFLOPS, LaplaceBoundGOPS, LaplaceBoundTime));
  }

  return (Timing());
} // PrintGFLOPSLaplaceBound

}    // namespace oclpdekernels
}
}

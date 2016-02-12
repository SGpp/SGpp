// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/pde/basis/common/OCLLaplaceInner.hpp>

using namespace SGPP::base;

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace parallel {
namespace oclpdekernels {

cl_kernel LaplaceInnerKernel[NUMDEVS];
cl_mem d_ptrLambdaInner[NUMDEVS];
double MultTimeLaplaceInner = 0.0;
double ReduTimeLaplaceInner = 0.0;
double CounterLaplaceInner = 0.0;
double LaplaceInnerStartupTime = 0.0;
double LaplaceInnerExecTime = 0.0;
double LaplaceInnerAllReduceTime = 0.0;
double LaplaceInnerExecStartTime = 0.0;
double LaplaceInnerExecMultTime = 0.0;
double LaplaceInnerExecReduceTime = 0.0;
double LaplaceInnerExecEndTime = 0.0;
SGppStopwatch* myStopwatch = new SGppStopwatch();
double* LaplaceInnerExecAll = new double[186];
double* LaplaceInnerProfiling  = new double[186];
double* LaplaceInnerWaiting = new double[186];


void SetLambdaBufferLaplaceInner(REAL* ptrLambda,
                                 size_t localdim) {

  for (size_t i = 0; i < 186; i++ ) {
    LaplaceInnerExecAll[i] = 0.0;
  }

  for (size_t i = 0; i < 186; i++ ) {
    LaplaceInnerProfiling[i] = 0.0;
  }

  for (size_t i = 0; i < 186; i++ ) {
    LaplaceInnerWaiting[i] = 0.0;
  }

  size_t lambda_size = localdim;
  cl_int ciErrNum = CL_SUCCESS;

  for (unsigned int i = 0; i < num_devices; ++i) {
    d_ptrLambdaInner[i] = clCreateBuffer(context,
                                         CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                         lambda_size * sizeof(REAL),
                                         ptrLambda, &ciErrNum);
    oclCheckErr(ciErrNum, "clCreateBuffer ptrLambda");
  }
}

std::string InnerGradientFunction() {
  std::stringstream stream_program_src;
  stream_program_src << "REAL gradient(REAL i_level,				"  << std::endl;
  stream_program_src << "      REAL i_index, 					"  << std::endl;
  stream_program_src << "      REAL j_level, 					"  << std::endl;
  stream_program_src << "      REAL j_index,					"  << std::endl;
  stream_program_src << "      REAL lcl_q_inv)					"  << std::endl;
  stream_program_src << "{							"  << std::endl;
  stream_program_src << "REAL grad;						"  << std::endl;
  stream_program_src <<
                     "ulong doGrad = (ulong)((i_level == j_level) && (i_index == j_index));	"  <<
                     std::endl;
  stream_program_src << "grad = select(0.0, i_level * 2.0 * lcl_q_inv,doGrad);			"
                     << std::endl;
  stream_program_src << "return grad;								"  << std::endl;
  stream_program_src << "}									"  << std::endl;
  return stream_program_src.str();
}

std::string LaplaceInnerHeader() {
  std::stringstream stream_program_src;
  stream_program_src << "__kernel void multKernel(__global  REAL* ptrLevel,			"
                     << std::endl;
  stream_program_src << "			 	__constant  REAL* ptrLevelIndexLevelintcon,	"  <<
                     std::endl;
  stream_program_src << "			 	__global  REAL* ptrIndex, 			"  << std::endl;
  stream_program_src << "			 	__global  REAL* ptrLevel_int, 			"  << std::endl;
  stream_program_src << "			 	__global  REAL* ptrAlpha, 			"  << std::endl;
  stream_program_src << "			 	__global  REAL* ptrParResult, 			"  << std::endl;
  stream_program_src << "			 	__constant  REAL* ptrLcl_q,			"  << std::endl;
  stream_program_src << "			 	__constant  REAL* ptrLambda,			"  << std::endl;
  stream_program_src << "			 	ulong overallMultOffset,   			"  << std::endl;
  stream_program_src << "                        	ulong ConstantMemoryOffset) 			"
                     << std::endl;
  stream_program_src << "{									"  << std::endl;
  return stream_program_src.str();
}


void CompileLaplaceInner(int id, std::string kernel_src, cl_kernel* kernel) {
  cl_int err = CL_SUCCESS;

  std::string source2 = LaplaceInnerHeader();
  std::string l2dotfunction = InnerLTwoDotFunction();
  std::string gradientfunction = InnerGradientFunction();
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
  stream_program_src << "__local REAL alphaTemp[" << LSIZE << "];" << std::endl;
  cl_uint sizec;
  clGetDeviceInfo(device_ids[0], CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV,
                  sizeof(cl_uint), &sizec, 0);


  size_t keplerUseRegister = (sizec >= 3);

  size_t fermiUseRegister = (dims <= 5);

#ifdef USEOCL_INTEL
  keplerUseRegister = 1;
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
                     "alphaTemp[get_local_id(0)]   = ptrAlpha[get_global_id(0) + overallMultOffset + ConstantMemoryOffset*LSIZE];"
                     << std::endl;

  stream_program_src << "REAL res = 0.0;" << std::endl;
  const char* typeREAL = "";

  if (keplerUseRegister) {
    const char* typeREAL2 = "REAL ";

    for (size_t d_inner = 0; d_inner < dims; d_inner++) {
      stream_program_src << typeREAL2 << "i_level_" << d_inner <<
                         " =         ptrLevel[" << d_inner* storageSizePadded <<
                         " + (get_global_id(1))*LSIZE];" << std::endl;
      stream_program_src << typeREAL2 << "i_index_" << d_inner <<
                         " =         ptrIndex[" << d_inner* storageSizePadded <<
                         " + (get_global_id(1))*LSIZE];" << std::endl;
      stream_program_src << typeREAL2 << "i_level_int_" << d_inner <<
                         " =     ptrLevel_int[" << d_inner* storageSizePadded <<
                         " + (get_global_id(1))*LSIZE];" << std::endl;
    }
  }


  stream_program_src << "for (unsigned k = 0; k < " << LSIZE << "; k++) {" <<
                     std::endl;

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
      stream_program_src << typeREAL2 << "i_level" << registerIdStr <<
                         " =         ptrLevel[" << d_inner* storageSizePadded <<
                         " + (get_global_id(1))*LSIZE];" << std::endl;
      stream_program_src << typeREAL2 << "i_index" << registerIdStr <<
                         " =         ptrIndex[" << d_inner* storageSizePadded <<
                         " + (get_global_id(1))*LSIZE];" << std::endl;
      stream_program_src << typeREAL2 << "i_level_int" << registerIdStr <<
                         " =     ptrLevel_int[" << d_inner* storageSizePadded <<
                         " + (get_global_id(1))*LSIZE];" << std::endl;
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
                         " + get_local_id(0)] = (l2dot(i_level" << registerIdStr << ",\
             i_index" << registerIdStr << ",\
             i_level_int" << registerIdStr << ", \
            j_levelreg,j_indexreg,j_level_intreg," << "ptrLcl_q[" <<
                         (d_inner + 1) * 2 - 2 << "]));" << std::endl;

      stream_program_src << "gradTemp[" << d_inner* LSIZE <<
                         " + get_local_id(0)] = (gradient(i_level" << registerIdStr << ",\
               i_index" << registerIdStr << ", \
               j_levelreg,j_indexreg," << "ptrLcl_q[" <<
                         (d_inner + 1) * 2 - 1 << "]));" << std::endl;
    } else {
      stream_program_src << "REAL l2dotreg_" << d_inner << " = \
	                 (l2dot(i_level" << registerIdStr << ",\
	                        i_index" << registerIdStr << ",\
	                        i_level_int" << registerIdStr << ",\
	                        j_levelreg,j_indexreg,j_level_intreg," << "ptrLcl_q["
                         << (d_inner + 1) * 2 - 2 << "]));" << std::endl;

      stream_program_src << "REAL gradreg_" << d_inner << " = \
	               (gradient(i_level" << registerIdStr << ",\
	                         i_index" << registerIdStr << ", \
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
                     "ptrParResult[((get_group_id(0)+ConstantMemoryOffset)*get_global_size(1) + get_global_id(1))*"
                     << LSIZE << " + get_local_id(0)] = res; " << std::endl;
  stream_program_src << "}" << std::endl;

  std::string program_src = stream_program_src.str();
  const char* source3 = program_src.c_str();
#if PRINTOCL
  std::cout <<  source3 << std::endl;
#endif
  cl_program program = clCreateProgramWithSource(context, 1,
                       (const char**)&source3, NULL, &err);
  oclCheckErr(err, "clCreateProgramWithSource");
  char buildOptions[256];
  int ierr = snprintf(buildOptions, sizeof(buildOptions),
                      "-DSTORAGE=%lu -DSTORAGEPAD=%lu -DNUMGROUPS=%lu -DDIMS=%lu -cl-finite-math-only -cl-fast-relaxed-math ",
                      storageSize, storageSizePadded, num_groups, dims);

  if (ierr < 0) {
    printf("Error in Build Options");
    exit(-1);
  }

  err = clBuildProgram(program, 0, NULL, buildOptions, NULL, NULL);

  if (err != CL_SUCCESS) {
    std::cout <<
              "OCL Error: compileMultKernel2DIndexesJInReg2 OpenCL Build Error. Error Code: "
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

}


void CompileLaplaceInnerKernels() {
  for (unsigned int i = 0; i < num_devices; ++i) {
    CompileLaplaceInner(i, "multKernel", LaplaceInnerKernel);
    CompileReduceInner(i, "ReduceInnerKernel", ReduceInnerKernel);
  }
}

void SetArgumentsLaplaceInner() {
  cl_int ciErrNum = CL_SUCCESS;
  int counter = 0;

  for (unsigned int i = 0; i < num_devices; ++i) {
    ciErrNum |= clSetKernelArg(LaplaceInnerKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrLevelInner[i]);
    ciErrNum |= clSetKernelArg(LaplaceInnerKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrLevelIndexLevelintconInner[i]);

    ciErrNum |= clSetKernelArg(LaplaceInnerKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrIndexInner[i]);
    ciErrNum |= clSetKernelArg(LaplaceInnerKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrLevel_intInner[i]);

    ciErrNum |= clSetKernelArg(LaplaceInnerKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrAlphaInner[i]);
    ciErrNum |= clSetKernelArg(LaplaceInnerKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrParResultInner[i]);
    ciErrNum |= clSetKernelArg(LaplaceInnerKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrLcl_qInner[i]);
    ciErrNum |= clSetKernelArg(LaplaceInnerKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrLambdaInner[i]);
    oclCheckErr(ciErrNum, "clSetKernelArg1 Kernel Construct");

    counter = 0;
    ciErrNum |= clSetKernelArg(ReduceInnerKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrResultInner[i]);
    ciErrNum |= clSetKernelArg(ReduceInnerKernel[i], counter++, sizeof(cl_mem),
                               (void*) &d_ptrParResultInner[i]);
    counter = 0;
    oclCheckErr(ciErrNum, "clSetKernelArg1 Kernel Construct");
  }
}


void CleanUpLaplaceInner() {
  if (!isFirstTimeLaplaceInner) {
    cl_int ciErrNum = CL_SUCCESS;

    for (unsigned int i = 0; i < num_devices; i++)  {
      ciErrNum |= clReleaseMemObject(d_ptrLambdaInner[i]);
      ciErrNum |= clReleaseKernel( LaplaceInnerKernel[i] );
      oclCheckErr(ciErrNum, "clReleaseMemObject & clReleaseKernel");
    }
  }
}


Timing PrintGFLOPSLaplaceInner() {

  double StorD = (double)storageSize;
  double dimsD = (double)dims;
  double LaplaceInnerTime = (MultTimeLaplaceInner + ReduTimeLaplaceInner) *
                            (1e-9);
  LaplaceInnerTime = LaplaceInnerExecTime;

  if (!isFirstTimeLaplaceInner) {
    std::cout << "Overall Runtime of AllReduceLaplaceInner: " <<
              LaplaceInnerAllReduceTime << std::endl;
    std::cout << "Time for GPU-calculation of LaplaceInner: " << LaplaceInnerTime <<
              std::endl;

    //    for(size_t i = 0; i < 186; i++) {
    //      std::cout << "Time for GPU-calculation of LaplaceInner["<<i<<"]: " << LaplaceInnerExecAll[i] << std::endl;

    //    }

  }

  // Laplace Inner
  double numberOfFlop = StorD * StorD * (dimsD * (26.0 + dimsD));

  if (!isFirstTimeLaplaceInner) {
    double LaplaceInnerGFLOPS = ((numberOfFlop * CounterLaplaceInner) /
                                 LaplaceInnerTime) * (1e-9);
    double LaplaceInnerOPS = StorD * StorD * (dimsD * (16.0));
    double LaplaceInnerGOPS = ((LaplaceInnerOPS * CounterLaplaceInner) /
                               LaplaceInnerTime) * (1e-9);
    return (Timing(LaplaceInnerGFLOPS, LaplaceInnerGOPS, LaplaceInnerTime));
  }

  return (Timing());
}

}
}
}


// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/pde/basis/common/OCLPDEInner.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
using namespace SGPP::base;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    namespace oclpdekernels {

      extern cl_kernel LaplaceInnerKernel[NUMDEVS];
      extern cl_mem d_ptrLambdaInner[NUMDEVS];
      extern double MultTimeLaplaceInner;
      extern double ReduTimeLaplaceInner;
      extern double CounterLaplaceInner;
      // Timing
      extern double LaplaceInnerStartupTime;
      extern double LaplaceInnerExecTime;
      extern double LaplaceInnerAllReduceTime;
      extern double LaplaceInnerExecStartTime;
      extern double LaplaceInnerExecMultTime;
      extern double LaplaceInnerExecReduceTime;
      extern double LaplaceInnerExecEndTime;

      extern SGppStopwatch* myStopwatch;
      extern double* LaplaceInnerExecAll;
      extern double* LaplaceInnerProfiling;
      extern double* LaplaceInnerWaiting;

      /// Allocates extra buffer for the Lambda parameter needed for the Laplace Operator.
      void SetLambdaBufferLaplaceInner(REAL* ptrLambda,
                                       size_t localdim);
      /// Returns the string with the OpenCL code for the function declaration for the Laplace operator on the inner grid.
      std::string LaplaceInnerHeader();
      /// Returns the string with the OpenCL code for the function for the Gradient function on the inner grid.
      std::string InnerGradientFunction();

      /// Generates and compiles the OpenCL code for the function for the Laplace operator on the inner grid.
      void CompileLaplaceInner(int id, std::string kernel_src, cl_kernel* kernel);


      /// Compiles all kernels pertaining to the Laplace operator (Laplace kernel, Reduction kernel) on inner grids.
      void CompileLaplaceInnerKernels();

      /// Sets arguments for all kernels pertaining to the Laplace operator on inner grids.
      void SetArgumentsLaplaceInner();

    }
  }
}


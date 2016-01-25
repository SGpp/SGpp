// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OCLPDEBound.hpp"
#include <sgpp/base/tools/SGppStopwatch.hpp>
using namespace SGPP::base;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    namespace oclpdekernels {

      extern cl_kernel LTwoDotBoundKernel[NUMDEVS];
      extern double MultTimeLTwoDotBound;
      extern double ReduTimeLTwoDotBound;
      extern double CounterLTwoDotBound;
      extern double LTwoDotBoundStartupTime;
      extern double LTwoDotBoundExecTime;
      extern double LTwoDotBoundAllReduceTime;
      extern SGppStopwatch* myStopwatch;

      /// Generates and compiles the OpenCL code for the function for the LTwoDot operator on the boundary grid.
      void CompileLTwoDotBound(int id, std::string kernel_src, cl_kernel* kernel);

      /// Compiles all kernels pertaining to the LTwoDot operator (LTwoDot kernel, Reduction kernel) on boundary grids.
      void CompileLTwoDotBoundKernels();

      /// Sets arguments for all kernels pertaining to the LTwoDot operator on boundary grids.
      void SetArgumentsLTwoDotBound();


    } // namespace oclpdekernels
  }
}

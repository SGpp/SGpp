#include "OCLLaplaceInner.hpp"

namespace sg {
  namespace parallel {
    namespace oclpdekernels {

      extern cl_kernel LTwoDotLaplaceInnerKernel[NUMDEVS];


      extern double LTwoDotLaplaceInnerAllReduceTime;
      extern double LTwoDotLaplaceInnerExecTime;
      extern double LTwoDotLaplaceInnerReduceTime;
      extern double LTwoDotLaplaceInnerExecTimeFirst;
      extern double LTwoDotLaplaceInnerExecTimeLast;
      extern double CounterLTwoDotLaplaceInner;
      extern double CounterLTwoDotLaplaceInnerFirst;
      extern double CounterLTwoDotLaplaceInnerRest;
      extern double CounterLTwoDotLaplaceInnerLast;
      extern double LTwoDotLaplaceInnerStartupTime;
      extern double * LTwoDotLaplaceInnerAll;
      extern double * LTwoDotLaplaceInnerProfiling;
      extern double * LTwoDotLaplaceInnerWaiting;
      extern double LTwoDotLaplaceInnerProfilingAcc ;
      extern double LTwoDotLaplaceInnerProfilingWait;

      /// Returns the string with the OpenCL code for the function declaration for the combined LTwoDot + Laplace operator on the inner grid.
      std::string LTwoDotLaplaceInnerHeader();

      /// Generates and compiles the OpenCL code for the function for the combined LTwoDot+Laplace operator on the inner grid.
      void CompileLTwoDotLaplace(int id, std::string kernel_src, cl_kernel* kernel);

      /// Compiles all kernels pertaining to the combined LTwoDot + Laplace operator on inner grids.
      void CompileLTwoDotLaplaceInnerKernels();

      /// Sets arguments for all kernels pertaining to the combined LTwoDot+Laplace operator on inner grids.
      void SetArgumentsLTwoDotLaplaceInner();

      /// Deallocates all data pertaining to the combined LTwoDot+Laplace Operator working on the inner grid  
      void CleanUpLTwoDotLaplaceInner();
    }
  }
}


#include "OCLLaplaceBound.hpp"

namespace sg {
  namespace parallel {
    namespace oclpdekernels {

      extern cl_kernel LTwoDotLaplaceBoundKernel[NUMDEVS];
      extern REAL TimestepCoeff;


      extern double LTwoDotLaplaceBoundAllReduceTime;
      extern double LTwoDotLaplaceBoundExecTime;
      extern double CounterLTwoDotLaplaceBound;
      extern double LTwoDotLaplaceBoundStartupTime;
      extern double* LTwoDotLaplaceBoundAll;

      /// Returns the string with the OpenCL code for the function declaration for the combined LTwoDot + Laplace operator on the boundary grid.
      std::string LTwoDotLaplaceBoundHeader();

      /// Generates and compiles the OpenCL code for the function for the combined LTwoDot+Laplace operator on the boundary grid.
      void CompileLTwoDotLaplaceBound(int id, std::string kernel_src, cl_kernel* kernel);

      /// Compiles all kernels pertaining to the combined LTwoDot + Laplace operator on boundary grids.
      void CompileLTwoDotLaplaceBoundKernels();

      /// Sets arguments for all kernels pertaining to the combined LTwoDot+Laplace operator on boundary grids.
      void SetArgumentsLTwoDotLaplaceBound();


    }
  }
}


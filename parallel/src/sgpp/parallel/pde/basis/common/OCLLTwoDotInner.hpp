#include "OCLPDEInner.hpp"
#include <sgpp/base/tools/SGppStopwatch.hpp>
using namespace SGPP::base;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    namespace oclpdekernels {


      extern cl_kernel LTwoDotInnerKernel[NUMDEVS];
      extern double MultTimeLTwoDotInner;
      extern double ReduTimeLTwoDotInner;
      extern double CounterLTwoDotInner;
      extern double LTwoDotInnerStartupTime;
      extern double LTwoDotInnerExecTime;
      extern double LTwoDotInnerAllReduceTime;
      extern SGppStopwatch* myStopwatch;

      /// Generates and compiles the OpenCL code for the function for the LTwoDot operator on the inner grid.
      void compileLTwoDotInner(int id, std::string kernel_src, cl_kernel* kernel);

      /// Compiles all kernels pertaining to the LTwoDot operator (LTwoDot kernel, Reduction kernel) on inner grids.
      void CompileLTwoDotInnerKernels();

      /// Sets arguments for all kernels pertaining to the LTwoDot operator on inner grids.
      void SetArgumentsLTwoDotInner();

      /// Deallocates all data pertaining to LTwoDot Operator working on the inner grid
      void CleanUpLTwoDotInner();

    } // namespace parallel
  }
}

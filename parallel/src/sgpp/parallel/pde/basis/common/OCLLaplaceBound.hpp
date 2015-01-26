#include "OCLPDEBound.hpp"
#include <sgpp/base/tools/SGppStopwatch.hpp>

using namespace SGPP::base;
#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    namespace oclpdekernels {

      extern cl_kernel LaplaceBoundKernel[NUMDEVS];
      extern cl_mem d_ptrLambdaBound[NUMDEVS];
      extern double MultTimeLaplaceBound;
      extern double ReduTimeLaplaceBound;
      extern double CounterLaplaceBound;
      extern double LaplaceBoundStartupTime;
      extern double LaplaceBoundExecTime;
      extern double LaplaceBoundAllReduceTime;
      extern SGppStopwatch* myStopwatch;

      /// Allocates extra buffer for the Lambda parameter needed for the Laplace Operator.
      void SetLambdaBufferLaplaceBound(REAL* ptrLambda,
                                       size_t localdim);
      /// Returns the string with the OpenCL code for the function declaration for the Laplace operator on the boundary grid.
      std::string LaplaceBoundHeader();
      /// Returns the string with the OpenCL code for the function for the Gradient function on the boundary grid.
      std::string BoundGradientFunction();

      /// Generates and compiles the OpenCL code for the function for the Laplace operator on the boundary grid.
      void CompileLaplaceBound(int id, std::string kernel_src, cl_kernel* kernel);

      /// Compiles all kernels pertaining to the Laplace operator (Laplace kernel, Reduction kernel) on boundary grids.
      void CompileLaplaceBoundKernels();

      /// Sets arguments for all kernels pertaining to the Laplace operator on boundary grids.
      void SetArgumentsLaplaceBound();


    }    // namespace oclpdekernels
  }
}

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/pde/basis/common/OCLPDEBound.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
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
extern sgpp::base::SGppStopwatch* myStopwatch;

/// Allocates extra buffer for the Lambda parameter needed for the Laplace Operator.
void SetLambdaBufferLaplaceBound(REAL* ptrLambda, size_t localdim);
/// Returns the string with the OpenCL code for the function declaration for the Laplace operator on
/// the boundary grid.
std::string LaplaceBoundHeader();
/// Returns the string with the OpenCL code for the function for the Gradient function on the
/// boundary grid.
std::string BoundGradientFunction();

/// Generates and compiles the OpenCL code for the function for the Laplace operator on the boundary
/// grid.
void CompileLaplaceBound(int id, std::string kernel_src, cl_kernel* kernel);

/// Compiles all kernels pertaining to the Laplace operator (Laplace kernel, Reduction kernel) on
/// boundary grids.
void CompileLaplaceBoundKernels();

/// Sets arguments for all kernels pertaining to the Laplace operator on boundary grids.
void SetArgumentsLaplaceBound();

}  // namespace oclpdekernels
}  // namespace parallel
}  // namespace sgpp

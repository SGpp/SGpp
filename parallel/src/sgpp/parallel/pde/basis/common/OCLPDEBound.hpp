// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_MPI
#ifdef USEOCL
#include <mpi.h>
#endif
#endif
#include <sgpp/parallel/pde/basis/common/OCLPDEKernels.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace parallel {
namespace oclpdekernels {

extern cl_kernel ReduceBoundKernel[NUMDEVS];
extern cl_mem d_ptrLevelBound[NUMDEVS];
extern cl_mem d_ptrIndexBound[NUMDEVS];
extern cl_mem d_ptrLevel_intBound[NUMDEVS];
extern cl_mem d_ptrResultBound[NUMDEVS];
extern cl_mem d_ptrResultPinnedBound[NUMDEVS];

extern cl_mem d_ptrParResultBound[NUMDEVS];
extern cl_mem d_ptrAlphaBound[NUMDEVS];
extern cl_mem d_ptrLevelIndexLevelintconBound[NUMDEVS];  // constant memory buffer holding all three
                                                         // components
extern cl_mem d_ptrLcl_qBound[NUMDEVS];                  // Also holds q_inverse

extern unsigned* offsetBound;
extern REAL* ptrResultBound;

extern REAL* ptrLevelTBound;
extern REAL* ptrIndexTBound;
extern REAL* ptrLevel_intTBound;
extern REAL* ptrParResultBound;
extern REAL* ptrAlphaEndBound;
extern REAL*
    ptrLevelIndexLevelintBound;  // for the constant memory buffer holding all three components
extern REAL* ptrLcl_qBound;      // Also holds q_inverse
extern REAL* ptrResultTempBound;
extern REAL* ptrResultZeroBound;
extern REAL* ptrResultPinnedBound;

extern size_t storageSizeBound;
extern size_t storageSizePaddedBound;
extern size_t num_groupsBound;
extern size_t max_buffer_sizeBound;
extern size_t par_result_sizeBound;
extern size_t lcl_q_sizeBound;
extern size_t alphaend_sizeBound;

extern size_t storageInnerSizePaddedBound;
extern size_t storageInnerSizeBound;
extern size_t InnerSizeBound;
extern size_t Inner_num_groupsBound;
extern size_t Inner_par_result_sizeBound;
extern size_t Inner_par_result_max_sizeBound;
extern size_t Inner_result_sizeBound;
extern size_t par_result_max_sizeBound;

extern size_t constant_mem_size;
extern size_t constant_buffer_size;
extern size_t constant_buffer_iterations;

extern size_t isFirstTimeLaplaceBound;
extern size_t isFirstTimeLTwoDotBound;
extern size_t isFirstTimeLTwoDotLaplaceBound;

/// Returns the string with the OpenCL code for the Reduction kernel for the operators on the
/// boundary grid.
std::string ReduceBoundKernelStr();
/// Compiles the OpenCL code for the Reduction kernel for the operators on the boundary grid and
/// saves it in kernel[id]. kernel_src must match the name of the OpenCL function.
void CompileReduceBound(int id, std::string kernel_src, cl_kernel* kernel);
/// Returns the string with the OpenCL code for the LTwoDot function for the operators on the
/// boundary grid.
std::string BoundLTwoDotFunction();
/// Allocates and initializes the main part of the buffers needed by the OpenCL code for the
/// operators on the boundary grid.
void SetBuffersBound(REAL* ptrLevel, REAL* ptrIndex, REAL* ptrLevel_int, size_t localStorageSize,
                     size_t localdim, sgpp::base::GridStorage* storage);
#ifdef USE_MPI
extern int* MPIOffsetListBound;
extern int* MPISizeListBound;
/// Creates partitions of the work load along the j-direction for the operators on the boundary
/// grid.
void SetUpMPIBound();
/// Computes the sum of all subresults result on all processes. Then it distributes results for the
/// inner grid back into the complete grid containing the boundary
void MPI_ShareResultAllReduceBound(sgpp::base::DataVector& result);
#endif
/// Deallocates all data pertaining to the Laplace operator on the boundary grid
void CleanUpLaplaceBound();
/// Deallocates all data pertaining to the LTwoDot operator on the boundary grid
void CleanUpLTwoDotBound();
/// Deallocates all data pertaining to the combined LTwoDot+Laplace Operator working on the boundary
/// grid
void CleanUpLTwoDotLaplaceBound();
}  // namespace oclpdekernels
}  // namespace parallel
}  // namespace sgpp

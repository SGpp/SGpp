// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OCLPDEKERNELS_HEADER
#define OCLPDEKERNELS_HEADER

#ifdef SEEK_SET
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#endif
#ifdef USE_MPI
#ifdef USEOCL
#include <mpi.h>
#endif
#endif

#ifdef USEOCL
#include <CL/cl.h>
#include <CL/cl_ext.h>
#endif
#include <string.h>
#include <malloc.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedGrid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <limits>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace parallel {

#define REAL double

/**
 * Implementation for
 *   - linear functions of Laplace Operation, linear grids without boundaries
 *   - linear functions of Laplace Operation, linear grids with boundaries
 *   - linear functions of LTwoDot Operation, linear grids without boundaries
 *   - linear functions of LTwoDot Operation, linear grids with boundaries
 *   - linear functions of the combined LTwoDot+Laplace Operation, linear grids without boundaries
 *   - linear functions of the combined LTwoDot+Laplace Operation, linear grids with boundaries
 *   - creation of the full matrix of the LTwoDot Operation, linear grids without boundaries
 * Internal data of the class is kept entirely in the global memory space within the SGPP::parallel::oclpdekernels namespace.
 */
class OCLPDEKernels {

 public:
  OCLPDEKernels() {
  }
  ;
  /**
   * Carries out the multiplication of a vector with the matrix for the Laplace operator on the inner grid.
   *
   * @param alpha SGPP::base::DataVector that contains the ansatzfunctions' coefficients
   * @param result SGPP::base::DataVector into which the result of the space discretization operation is stored
   * @param lcl_q Array with the interval width for each dimension
   * @param lcl_q_inv Array with the reciprocal interval width for each dimension
   * @param ptrLevel Array with the level of each ansatzfunction.
   * @param ptrIndex Array with the Index of each ansatzfunction.
   * @param ptrLevel_int Array with the level for the integral of each ansatzfunction.
   * @param ptrLambda the lambda parameter which is needed in some cases (Black-Scholes) to modify the dimensional local values.
   * @param argStorageSize Size of the outer dimension of the ptrLevel, ptrIndex and ptrLevel_int arrays.
   * @param argStorageDim Size of the inner dimension of the ptrLevel, ptrIndex and ptrLevel_int arrays.
   * @param storage Grid which is to be used for evaluation.
   */
  void RunOCLKernelLaplaceInner(SGPP::base::DataVector& alpha,
                                SGPP::base::DataVector& result,
                                REAL* lcl_q,
                                REAL* lcl_q_inv,
                                REAL* ptrLevel,
                                REAL* ptrIndex,
                                REAL* ptrLevel_int,
                                REAL* ptrLambda, size_t argStorageSize, size_t argStorageDim,
                                SGPP::base::GridStorage* storage);

  /**
   * Carries out the multiplication of a vector with the matrix for the Laplace operator on the boundary grid.
   *
   */
  void RunOCLKernelLaplaceBound(SGPP::base::DataVector& alpha,
                                SGPP::base::DataVector& result,
                                REAL* lcl_q,
                                REAL* lcl_q_inv,
                                REAL* ptrLevel,
                                REAL* ptrIndex,
                                REAL* ptrLevel_int,
                                REAL* ptrLambda, size_t argStorageSize, size_t argStorageDim,
                                SGPP::base::GridStorage* storage);

  /**
   * Carries out the multiplication of a vector with the matrix for the LTwoDot operator on the inner grid.
   *
   */
  void RunOCLKernelLTwoDotInner(SGPP::base::DataVector& alpha,
                                SGPP::base::DataVector& result,
                                REAL* lcl_q,
                                REAL* ptrLevel,
                                REAL* ptrIndex,
                                REAL* ptrLevel_int, size_t argStorageSize, size_t argStorageDim,
                                SGPP::base::GridStorage* storage);

  /**
   * Carries out the multiplication of a vector with the matrix for the LTwoDot operator on the boundary grid.
   *
   */
  void RunOCLKernelLTwoDotBound(SGPP::base::DataVector& alpha,
                                SGPP::base::DataVector& result,
                                REAL* lcl_q,
                                REAL* ptrLevel,
                                REAL* ptrIndex,
                                REAL* ptrLevel_int, size_t argStorageSize, size_t argStorageDim,
                                SGPP::base::GridStorage* storage);

  /**
   * Carries out the multiplication of a vector with the matrix for the combined LTwoDot+Laplace operator on the inner grid.
   *
   */
  void RunOCLKernelLTwoDotLaplaceInner(SGPP::base::DataVector& alpha,
                                       SGPP::base::DataVector& result,
                                       REAL* lcl_q,
                                       REAL* lcl_q_inv,
                                       REAL* ptrLevel,
                                       REAL* ptrIndex,
                                       REAL* ptrLevel_int,
                                       REAL* ptrLambda, size_t argStorageSize, size_t argStorageDim,
                                       SGPP::base::GridStorage* storage,
                                       REAL tsCoeff);

  void RunOCLKernelSymmetricLTwoDotLaplaceInner(SGPP::base::DataVector& alpha,
      SGPP::base::DataVector& result,
      REAL* lcl_q,
      REAL* lcl_q_inv,
      REAL* ptrLevel,
      REAL* ptrIndex,
      REAL* ptrLevel_int,
      REAL* ptrLambda, size_t argStorageSize, size_t argStorageDim,
      SGPP::base::GridStorage* storage,
      REAL tsCoeff);

  void RunOCLKernelLTwoDotLaplaceInnerCPU(SGPP::base::DataVector& alpha,
                                          SGPP::base::DataVector& result,
                                          REAL* lcl_q,
                                          REAL* lcl_q_inv,
                                          REAL* ptrLevel,
                                          REAL* ptrIndex,
                                          REAL* ptrLevel_int,
                                          REAL* ptrLambda, size_t argStorageSize, size_t argStorageDim,
                                          SGPP::base::GridStorage* storage,
                                          REAL tsCoeff);
  /**
   * Carries out the multiplication of a vector with the matrix for the combined LTwoDot+Laplace operator on the boundary grid.
   *
   */
  void RunOCLKernelLTwoDotLaplaceBound(SGPP::base::DataVector& alpha,
                                       SGPP::base::DataVector& result,
                                       REAL* lcl_q,
                                       REAL* lcl_q_inv,
                                       REAL* ptrLevel,
                                       REAL* ptrIndex,
                                       REAL* ptrLevel_int,
                                       REAL* ptrLambda, size_t argStorageSize, size_t argStorageDim,
                                       SGPP::base::GridStorage* storage,
                                       REAL tsCoeff);

  /**
   * Generates the full matrix for the combined LTwoDot operator on the inner grid.
   *
   */
  void RunOCLKernelGenAInner(SGPP::base::DataVector& alpha,
                             REAL* ptrA,
                             REAL* lcl_q,
                             REAL* lcl_q_inv,
                             REAL* ptrLevel,
                             REAL* ptrIndex,
                             REAL* ptrLevel_int,
                             REAL* ptrLambda, size_t argStorageSize, size_t argStorageDim,
                             SGPP::base::GridStorage* storage);

  /**
   * Deallocates data pertaining to all the operators.
   */
  void CleanUpGPU();
};
// Class OCLPDEKernels

namespace oclpdekernels {

#define NUMDEVS 1

#define MAXBYTES 2100000000
#ifdef USEOCL_MIC
#define LSIZE 16
#elif USEOCL_CPU
#define LSIZE 16
#else
#define LSIZE 64
#endif
#define MOD 186

#define TOTALTIMING 1
#define PROFILING 0
#define QUEUEPROF 0
#define PRINTOCL  0
#define PRINTBUFFERSIZES 0

extern cl_uint num_devices;
extern cl_uint num_platforms;
extern cl_platform_id platform_id;
extern cl_platform_id* platform_ids;
extern cl_device_id* device_ids;
extern cl_context context;
extern cl_command_queue command_queue[NUMDEVS];

extern size_t isVeryFirstTime;
extern size_t isCleanedUp;
extern double AverageGFLOPS;
extern double GPUNVDMAXFLOPS;
extern double GPUNVDMAXFLOPSHALF;
extern size_t padding_size;
extern size_t dims;
extern REAL TimestepCoeff;
extern size_t LastRunBound;

/// Write an error and stops program if an error occurs.
void oclCheckErr(cl_int err, const char* function);

/// Initializes the platform, CPU or GPU, and the command queues.
void StartUpGPU();

/// Class for passing timing results for individual Operators
class Timing {
 public:
  double GFLOPS;
  double GOPS;
  double time;
  Timing() :
    GFLOPS(0.0), GOPS(0.0), time(0.0) {
  }
  Timing(double iGFLOPS, double iGOPS, double itime) :
    GFLOPS(iGFLOPS), GOPS(iGOPS), time(itime) {
  }
};

/// Deallocates all data pertaining to operators working on the inner grid
void CleanUpInner();

/// Deallocates all data pertaining to operators working on the boundary grid
void CleanUpBound();

/// Prints and returns timing results for the Laplace Operator on the inner grid
Timing PrintGFLOPSLaplaceInner();
/// Prints and returns timing results for the Laplace Operator on the boundary grid
Timing PrintGFLOPSLaplaceBound();
/// Prints and returns timing results for the Laplace Operator on the inner grid
Timing PrintGFLOPSLTwoDotInner();
/// Prints and returns timing results for the Laplace Operator on the boundary grid
Timing PrintGFLOPSLTwoDotBound();
/// Prints and returns timing results for the Laplace Operator on the inner grid
Timing PrintGFLOPSLTwoDotLaplaceInner();
/// Prints and returns timing results for the Laplace Operator on the boundary grid
Timing PrintGFLOPSLTwoDotLaplaceBound();
/// returns the timing (nanoseconds) of a kernel execution given by GPUExecution
double AccumulateTiming(cl_event* GPUExecution, size_t gpuid);
/// returns the timing (nanoseconds) of a kernel waiting in the queue given by GPUExecution
double AccumulateWaiting(cl_event* GPUExecution, size_t gpuid);
}
}
}
#endif //OCLPDEKERNELS_HEADER

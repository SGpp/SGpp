// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/pde/basis/common/OCLLTwoDotLaplaceBound.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace parallel {
namespace oclpdekernels {

void ExecLTwoDotLaplaceBound(REAL* ptrAlpha, REAL* ptrResult, REAL* lcl_q, REAL* lcl_q_inv,
                             REAL* ptrLevel, REAL* ptrIndex, REAL* ptrLevel_int,
                             size_t argStorageSize, size_t argStorageDim, size_t MPIOffset,
                             size_t MPIglobalsize) {
  cl_int ciErrNum = CL_SUCCESS;
  cl_event GPUDone[NUMDEVS];
  cl_event GPUDoneLcl[NUMDEVS];
  cl_event GPUExecution[NUMDEVS];

  size_t idx = 0;

  for (size_t d_outer = 0; d_outer < dims; d_outer++) {
    ptrLcl_qBound[idx++] = lcl_q[d_outer];
    ptrLcl_qBound[idx++] = lcl_q_inv[d_outer];
  }

  for (size_t i = 0; i < num_devices; i++) {
    ciErrNum |= clEnqueueWriteBuffer(command_queue[i], d_ptrAlphaBound[i], CL_FALSE, 0,
                                     storageSizeBound * sizeof(REAL), ptrAlpha, 0, 0, &GPUDone[i]);
    ciErrNum |=
        clEnqueueWriteBuffer(command_queue[i], d_ptrLcl_qBound[i], CL_FALSE, 0,
                             lcl_q_sizeBound * sizeof(REAL), ptrLcl_qBound, 0, 0, &GPUDoneLcl[i]);
  }

  oclCheckErr(ciErrNum, "clEnqueueWriteBuffer L2225");
  clWaitForEvents(num_devices, GPUDone);
  clWaitForEvents(num_devices, GPUDoneLcl);

  for (size_t i = 0; i < num_devices; i++) {
    ciErrNum |= clEnqueueWriteBuffer(command_queue[i], d_ptrResultBound[i], CL_FALSE, 0,
                                     Inner_result_sizeBound * sizeof(REAL), ptrResultZeroBound, 0,
                                     0, &GPUDone[i]);
    ciErrNum |= clEnqueueWriteBuffer(
        command_queue[i], d_ptrAlphaBound[i], CL_FALSE, storageSizeBound * sizeof(REAL),
        alphaend_sizeBound * sizeof(REAL), ptrAlphaEndBound, 0, 0, &GPUExecution[i]);
  }

  oclCheckErr(ciErrNum, "clEnqueueWriteBuffer mult");
  clWaitForEvents(num_devices, GPUDone);

  {
    size_t multglobalworksize[2];
    size_t MPIglobal = MPIglobalsize;
    /* Find minimum size of first work-grid dimension for the total computation.
     * It is limited by
     * 1. size of parResult buffer
     * 2. size of discretization grid
     */
    size_t storageSizePaddedStep =
        std::min(storageSizePaddedBound / num_devices,
                 Inner_par_result_max_sizeBound / (storageSizePaddedBound)*LSIZE);
    multglobalworksize[0] = std::min(storageSizePaddedBound, storageSizePaddedStep);
    multglobalworksize[1] = storageInnerSizePaddedBound / LSIZE;

    // Maximum number of grid points in the first dimension of the discretization grid
    size_t multglobal = MPIglobal / num_devices;

    for (size_t overallMultOffset = 0; overallMultOffset < multglobal;
         overallMultOffset += std::min(multglobalworksize[0], multglobal - overallMultOffset)) {
      multglobalworksize[0] = std::min(multglobalworksize[0], multglobal - overallMultOffset);

      for (unsigned int i = 0; i < num_devices; i++) {
        /* offset in discretization grid based on three factors
         * 1. offset from MPI
         * 2. offset if parResult buffer is too large (only non-zero for large problems)
         * 3. offset from multi-GPU implementation, if used.
         */
        size_t overallMultOffset2 = MPIOffset + overallMultOffset + i * multglobal;
        ciErrNum |= clSetKernelArg(LTwoDotLaplaceBoundKernel[i], 8, sizeof(cl_ulong),
                                   (void*)&overallMultOffset2);
        oclCheckErr(ciErrNum, "clSetKernelArgL1660");
      }

      size_t constantglobalworksize[2];
      size_t constantlocalworksize[2];
      size_t constantglobal = multglobalworksize[0];
      /* Find minimum size of first work-grid dimension when we block for the constant memory.
       * It is limited by
       * 1. size of constant memory
       * 2. size of discretization grid
       * Usual value, e.g., 512.
       */
      constantglobalworksize[0] = std::min(constant_buffer_iterations, constantglobal);
      constantglobalworksize[1] = multglobalworksize[1];

      constantlocalworksize[0] = LSIZE;
      constantlocalworksize[1] = 1;

      for (size_t ConstantMemoryOffset = 0; ConstantMemoryOffset < constantglobal;
           ConstantMemoryOffset +=
           std::min(constantglobalworksize[0], constantglobal - ConstantMemoryOffset)) {
        constantglobalworksize[0] =
            std::min(constantglobalworksize[0], constantglobal - ConstantMemoryOffset);

        for (unsigned int i = 0; i < num_devices; i++) {
          // Write block of data needed by the kernel to constant memory.
          ciErrNum |= clEnqueueWriteBuffer(
              command_queue[i], d_ptrLevelIndexLevelintconBound[i], CL_FALSE, 0,
              constantglobalworksize[0] * 3 * dims * sizeof(REAL),
              ptrLevelIndexLevelintBound +
                  (MPIOffset + overallMultOffset + i * multglobal + ConstantMemoryOffset) * 3 *
                      dims,
              1, &GPUExecution[i], &GPUDone[i]);

          size_t jj = (ConstantMemoryOffset) / LSIZE;
          oclCheckErr(ciErrNum, "clEnqueueWriteBufferOCLLapBoundL2157");
          ciErrNum |= clSetKernelArg(LTwoDotLaplaceBoundKernel[i], 9, sizeof(cl_ulong), (void*)&jj);
          ciErrNum |=
              clSetKernelArg(LTwoDotLaplaceBoundKernel[i], 10, sizeof(REAL), (void*)&TimestepCoeff);

          oclCheckErr(ciErrNum, "clSetKernelArgOCLLTwoLapBoundL75");
          ciErrNum = clEnqueueNDRangeKernel(command_queue[i], LTwoDotLaplaceBoundKernel[i], 2, 0,
                                            constantglobalworksize, constantlocalworksize, 1,
                                            &GPUDone[i], &GPUExecution[i]);
          oclCheckErr(ciErrNum, "clEnqueueNDRangeKernel2195");
        }
      }

      // Perform reduction on the ParResult buffer
      size_t overallReduceOffset = 0;

      for (unsigned int i = 0; i < num_devices; i++) {
        // currently not used; always zero.
        ciErrNum |=
            clSetKernelArg(ReduceBoundKernel[i], 2, sizeof(cl_ulong), (void*)&overallReduceOffset);
        // Set the number of rows to sum in ParResult (or S)
        size_t newnum_groups = multglobalworksize[0] / LSIZE;
        ciErrNum |=
            clSetKernelArg(ReduceBoundKernel[i], 3, sizeof(cl_ulong), (void*)&newnum_groups);
        oclCheckErr(ciErrNum, "clSetKernelArgL1205");

        size_t reduceglobalworksize2[] = {multglobalworksize[1] * LSIZE, 1};
        size_t local2[] = {LSIZE, 1};
        ciErrNum |=
            clEnqueueNDRangeKernel(command_queue[i], ReduceBoundKernel[i], 2, 0,
                                   reduceglobalworksize2, local2, 0, NULL, &GPUExecution[i]);
        oclCheckErr(ciErrNum, "clEnqueueNDRangeKernel1213");
      }
    }
  }

  /* Read back the result to host, first to pinned buffer, then
   * to ptrResultBound which is used in the MPI_Allreduce function.
   */
  if (num_devices > 1) {
    for (unsigned int i = 0; i < num_devices; i++) {
      ciErrNum |= clEnqueueReadBuffer(
          command_queue[i], d_ptrResultBound[i], CL_FALSE, 0, storageInnerSizeBound * sizeof(REAL),
          ptrResultPinnedBound + i * storageInnerSizePaddedBound, 1, &GPUExecution[i], &GPUDone[i]);
    }

    oclCheckErr(ciErrNum, "clEnqueueReadBufferLapIL2145");
    ciErrNum |= clWaitForEvents(num_devices, GPUDone);
    oclCheckErr(ciErrNum, "clWaitForEventsLapIL2147");

    for (size_t i = 0; i < storageInnerSizeBound; i++) {
      ptrResultBound[i] = 0.0;
    }

    for (size_t j = 0; j < num_devices; j++) {
      for (size_t i = 0; i < storageInnerSizeBound; i++) {
        ptrResultBound[i] += ptrResultPinnedBound[j * storageInnerSizePaddedBound + i];
      }
    }
  } else {
    for (unsigned int i = 0; i < num_devices; i++) {
      ciErrNum |= clEnqueueReadBuffer(command_queue[i], d_ptrResultBound[i], CL_TRUE, 0,
                                      storageInnerSizeBound * sizeof(REAL), ptrResultPinnedBound, 1,
                                      &GPUExecution[i], &GPUDone[i]);
      oclCheckErr(ciErrNum, "clEnqueueReadBufferL161");
    }

    for (size_t i = 0; i < storageInnerSizeBound; i++) {
      ptrResultBound[i] = ptrResultPinnedBound[i];
    }
  }

#if TOTALTIMING
  CounterLTwoDotLaplaceBound += 1.0;
#endif

  for (size_t i = 0; i < num_devices; i++) {
    clReleaseEvent(GPUExecution[i]);
    clReleaseEvent(GPUDone[i]);
    clReleaseEvent(GPUDoneLcl[i]);
  }
}  // ExecLaplaceBound

}  // namespace oclpdekernels
using namespace oclpdekernels;
void OCLPDEKernels::RunOCLKernelLTwoDotLaplaceBound(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, REAL* lcl_q, REAL* lcl_q_inv,
    REAL* ptrLevel, REAL* ptrIndex, REAL* ptrLevel_int, REAL* ptrLambda, size_t argStorageSize,
    size_t argStorageDim, SGPP::base::GridStorage* storage, REAL tsCoeff) {
  TimestepCoeff = tsCoeff;
  myStopwatch->start();

  if (isVeryFirstTime) {
    StartUpGPU();
  }

  if (isFirstTimeLaplaceBound && isFirstTimeLTwoDotBound && isFirstTimeLTwoDotLaplaceBound) {
    SetBuffersBound(ptrLevel, ptrIndex, ptrLevel_int, argStorageSize, argStorageDim, storage);
    SetUpMPIBound();
  }

  if (isFirstTimeLaplaceBound && isFirstTimeLTwoDotLaplaceBound) {
    SetLambdaBufferLaplaceBound(ptrLambda, argStorageDim);
  }

  if (isFirstTimeLTwoDotLaplaceBound) {
    CompileLTwoDotLaplaceBoundKernels();
    SetArgumentsLTwoDotLaplaceBound();
    isVeryFirstTime = 0;
    isFirstTimeLTwoDotLaplaceBound = 0;
  }

  LTwoDotLaplaceBoundStartupTime += myStopwatch->stop();

  LastRunBound = 1;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (MPISizeListBound[myrank] != 0) {
    myStopwatch->start();
    ExecLTwoDotLaplaceBound(alpha.getPointer(), result.getPointer(), lcl_q, lcl_q_inv, ptrLevel,
                            ptrIndex, ptrLevel_int, argStorageSize, argStorageDim,
                            MPIOffsetListBound[myrank], MPISizeListBound[myrank]);
    double runtime = myStopwatch->stop();
    LTwoDotLaplaceBoundExecTime += runtime;

  } else {
    // Needed when the process has no work, otherwise it has false values.
    for (size_t i = 0; i < storageInnerSizeBound; i++) {
      ptrResultBound[i] = 0.0;
    }
  }

  myStopwatch->start();
  MPI_ShareResultAllReduceBound(result);
  LTwoDotLaplaceBoundAllReduceTime += myStopwatch->stop();

}  // RunOCLKernelLaplaceBound
}
}

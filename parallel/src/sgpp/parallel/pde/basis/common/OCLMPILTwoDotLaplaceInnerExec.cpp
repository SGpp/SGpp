#include "OCLLTwoDotLaplaceInner.hpp"

namespace sg {
  namespace parallel {
    namespace oclpdekernels {

      void ExecLTwoDotLaplaceInner(REAL* ptrAlpha,
                                   REAL* ptrResult,
                                   REAL* lcl_q,
                                   REAL* lcl_q_inv,
                                   REAL* ptrLevel,
                                   REAL* ptrIndex,
                                   REAL* ptrLevel_int,
                                   size_t argStorageSize,
                                   size_t argStorageDim,
                                   size_t MPIOffset,
                                   size_t MPIglobalsize) {

        cl_int ciErrNum = CL_SUCCESS;
        cl_event GPUDone[NUMDEVS];
        cl_event GPUDoneLcl[NUMDEVS];
        cl_event GPUExecution[NUMDEVS];

        size_t idx = 0;

        for (size_t d_outer = 0; d_outer < dims ; d_outer++) {
          ptrLcl_qInner[idx++] = lcl_q[d_outer];
          ptrLcl_qInner[idx++] = lcl_q_inv[d_outer];
        }

        for (size_t i = 0; i < storageSize; i++) {
          ptrAlphaPinnedInner[i] = ptrAlpha[i];
        }


        for (size_t i = 0; i < num_devices; i++) {

          ciErrNum |= clEnqueueWriteBuffer(command_queue[i], d_ptrAlphaInner[i], CL_FALSE, 0,
                                           storageSizePadded * sizeof(REAL), ptrAlphaPinnedInner, 0, 0, &GPUDone[i]);
          ciErrNum |= clEnqueueWriteBuffer(command_queue[i], d_ptrLcl_qInner[i], CL_FALSE, 0,
                                           lcl_q_size * sizeof(REAL), ptrLcl_qInner, 0, 0, &GPUDoneLcl[i]);
        }

        clWaitForEvents(num_devices, GPUDone);
        clWaitForEvents(num_devices, GPUDoneLcl);

        for (size_t i = 0; i < num_devices; i++) {
          ciErrNum |= clEnqueueWriteBuffer(command_queue[i], d_ptrResultInner[i],
                                           CL_FALSE, 0,
                                           storageSizePadded * sizeof(REAL),
                                           ptrResultZero, 0, 0, &GPUExecution[i]);
        }

        oclCheckErr(ciErrNum, "clEnqueueWriteBuffer mult");

        {
          size_t multglobalworksize[2];
          size_t MPIglobal = MPIglobalsize;
          /* Find minimum size of first work-grid dimension for the total computation.
           * It is limited by
           * 1. size of parResult buffer
           * 2. size of discretization grid
           */
          size_t storageSizePaddedStep = std::min(storageSizePadded  / num_devices, par_result_max_size / (storageSizePadded) * LSIZE);
          multglobalworksize[0] = std::min(storageSizePadded, storageSizePaddedStep);
          multglobalworksize[1] = storageSizePadded / LSIZE;
          // Maximum number of grid points in the first dimension of the discretization grid
          size_t multglobal = MPIglobal / num_devices;

          for (size_t overallMultOffset = 0; overallMultOffset < multglobal; overallMultOffset += std::min(multglobalworksize[0], multglobal - overallMultOffset)) {
            multglobalworksize[0] =  std::min(multglobalworksize[0], multglobal - overallMultOffset);

            for (unsigned int i = 0; i < num_devices; i++) {
              /* offset in discretization grid based on three factors
               * 1. offset from MPI
               * 2. offset if parResult buffer is too large (only non-zero for large problems)
               * 3. offset from multi-GPU implementation, if used.
               */
              size_t overallMultOffset2 = MPIOffset + overallMultOffset + i * multglobal;
              ciErrNum |= clSetKernelArg(LTwoDotLaplaceInnerKernel[i], 8, sizeof(cl_ulong), (void*) &overallMultOffset2);
            }

            oclCheckErr(ciErrNum, "clSetKernelArgL246");
            size_t constantglobalworksize[2];
            size_t constantlocalworksize[2];
            size_t constantglobal = multglobalworksize[0];
            /* Find minimum size of first work-grid dimension when we block for the constant memory.
             * It is limited by
             * 1. size of constant memory
             * 2. size of discretization grid
             * Usual value, e.g., 512.
             */
            constantglobalworksize[0] = std::min(constant_buffer_iterations_noboundary, constantglobal);
            constantglobalworksize[1] = multglobalworksize[1];
            constantlocalworksize[0] = LSIZE;
            constantlocalworksize[1] = 1;

            for (size_t ConstantMemoryOffset = 0; ConstantMemoryOffset < constantglobal; ConstantMemoryOffset += std::min(constantglobalworksize[0], constantglobal - ConstantMemoryOffset)) {
              constantglobalworksize[0] = std::min(constantglobalworksize[0], constantglobal - ConstantMemoryOffset);


              for (unsigned int i = 0; i < num_devices; i++) {
                // Write block of data needed by the kernel to constant memory.
                ciErrNum |= clEnqueueWriteBuffer(command_queue[i],
                                                 d_ptrLevelIndexLevelintconInner[i],
                                                 CL_FALSE, 0 ,
                                                 constantglobalworksize[0] * 3 * dims * sizeof(REAL),
                                                 ptrLevelIndexLevelintInner + (MPIOffset + overallMultOffset + i * multglobal + ConstantMemoryOffset) * 3 * dims,
                                                 1,
                                                 &GPUExecution[i] , &GPUDone[i]);
                oclCheckErr(ciErrNum, "clEnqueueWriteBufferL306");
                size_t jj = (ConstantMemoryOffset) / LSIZE;

                ciErrNum |= clSetKernelArg(LTwoDotLaplaceInnerKernel[i], 9, sizeof(cl_ulong), (void*) &jj);
                ciErrNum |= clSetKernelArg(LTwoDotLaplaceInnerKernel[i], 10, sizeof(REAL), (void*) &TimestepCoeff);

                oclCheckErr(ciErrNum, "clEnqueueWriteBufferL302");

                ciErrNum = clEnqueueNDRangeKernel(command_queue[i],
                                                  LTwoDotLaplaceInnerKernel[i],
                                                  2, 0,
                                                  constantglobalworksize,
                                                  constantlocalworksize,
                                                  1, &GPUDone[i], &GPUExecution[i]);
                oclCheckErr(ciErrNum, "clEnqueueNDRangeKernel2314");
              }
            }

            // Perform reduction on the ParResult buffer
            size_t overallReduceOffset = 0;//overallMultOffset*LSIZE;

            for (unsigned int i = 0; i < num_devices; i++) {
              // currently not used; always zero.
              ciErrNum |= clSetKernelArg(ReduceInnerKernel[i], 2, sizeof(cl_ulong), (void*) &overallReduceOffset);
              oclCheckErr(ciErrNum, "clSetKernelArgLapIL2199");
              // Set the number of rows to sum in ParResult (or S)
              size_t newnum_groups = multglobalworksize[0] / LSIZE ;
              ciErrNum |= clSetKernelArg(ReduceInnerKernel[i], 3, sizeof(cl_ulong), (void*) &newnum_groups);
              oclCheckErr(ciErrNum, "clSetKernelArgLapIL340");

              size_t reduceglobalworksize2[] = {multglobalworksize[1]* LSIZE, 1};
              size_t local2[] = {LSIZE, 1};
              ciErrNum |= clEnqueueNDRangeKernel(command_queue[i],
                                                 ReduceInnerKernel[i],
                                                 2, 0, reduceglobalworksize2,
                                                 local2,
                                                 0, NULL, &GPUExecution[i]);
              oclCheckErr(ciErrNum, "clEnqueueNDRangeKernel1213");
            }
          }
        }


        /* Read back the result to host, first to pinned buffer, then
         * to ptrResult which is used in the MPI_Allreduce function.
         */
        if (num_devices > 1) {

          for (unsigned int i = 0; i < num_devices; i++) {
            ciErrNum |= clEnqueueReadBuffer(command_queue[i],
                                            d_ptrResultInner[i], CL_FALSE, 0,
                                            storageSize * sizeof(REAL),
                                            ptrResultPinnedInner + i * storageSizePadded,
                                            1, &GPUExecution[i], &GPUDone[i]);

          }

          oclCheckErr(ciErrNum, "clEnqueueReadBufferLapIL2145");
          ciErrNum |= clWaitForEvents(num_devices, GPUDone);
          oclCheckErr(ciErrNum, "clWaitForEventsLapIL2147");

          for (size_t j = 0; j < num_devices; j++) {
            for (size_t i = 0; i < storageSize; i++) {
              ptrResult[i] += ptrResultPinnedInner[j * storageSizePadded + i];

            }
          }

        } else {
          for (unsigned int i = 0; i < num_devices; i++) {
            ciErrNum |= clEnqueueReadBuffer(command_queue[i],
                                            d_ptrResultInner[i],
                                            CL_TRUE, 0,
                                            storageSize * sizeof(REAL),
                                            ptrResultPinnedInner,
                                            1, &GPUExecution[i], &GPUDone[i]);
            oclCheckErr(ciErrNum, "clEnqueueReadBufferLapIL173");

          }

          for ( size_t i = 0; i < storageSize; i++) {
            ptrResult[i] = ptrResultPinnedInner[i];
          }
        }

#if TOTALTIMING
        CounterLTwoDotLaplaceInner += 1.0;
#endif

        for (size_t i = 0; i < num_devices; i++) {
          clReleaseEvent(GPUExecution[i]);
          clReleaseEvent(GPUDone[i]);
          clReleaseEvent(GPUDoneLcl[i]);
        }
      }

    }

    using namespace oclpdekernels;
    void OCLPDEKernels::RunOCLKernelLTwoDotLaplaceInner(sg::base::DataVector& alpha,
        sg::base::DataVector& result,
        REAL* lcl_q,
        REAL* lcl_q_inv,
        REAL* ptrLevel,
        REAL* ptrIndex,
        REAL* ptrLevel_int,
        REAL* ptrLambda,
        size_t argStorageSize,
        size_t argStorageDim,
        sg::base::GridStorage* storage,
        REAL tsCoeff) {

      TimestepCoeff = tsCoeff;
      myStopwatch->start();

      if (isVeryFirstTime) {
        StartUpGPU();
      }

      if (isFirstTimeLaplaceInner &&
          isFirstTimeLTwoDotInner &&
          isFirstTimeLTwoDotLaplaceInner) {
        SetBuffersInner(ptrLevel,
                        ptrIndex,
                        ptrLevel_int,
                        argStorageSize,
                        argStorageDim, storage);
        SetUpMPIInner();
      }

      if (isFirstTimeLaplaceInner && isFirstTimeLTwoDotLaplaceInner) {
        SetLambdaBufferLaplaceInner(ptrLambda,
                                    argStorageDim);

      }

      if (isFirstTimeLTwoDotLaplaceInner) {
        CompileLTwoDotLaplaceInnerKernels();
        SetArgumentsLTwoDotLaplaceInner();
        isVeryFirstTime = 0;
        isFirstTimeLTwoDotLaplaceInner = 0;
      }


      LTwoDotLaplaceInnerStartupTime += myStopwatch->stop();

      int myrank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

      if (MPISizeListInner[myrank] != 0) {
        myStopwatch->start();
        ExecLTwoDotLaplaceInner(alpha.getPointer(),
                                result.getPointer(),
                                lcl_q, lcl_q_inv,
                                ptrLevel, ptrIndex,
                                ptrLevel_int, argStorageSize,
                                argStorageDim,
                                MPIOffsetListInner[myrank],
                                MPISizeListInner[myrank]);
        double runtime = myStopwatch->stop();
        LTwoDotLaplaceInnerExecTime += runtime;
      }

      myStopwatch->start();
      MPI_CombineResultInner(result);
      LTwoDotLaplaceInnerAllReduceTime += myStopwatch->stop();
    }
  }
}


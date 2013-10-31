/* ****************************************************************************
* Copyright (C) 2012-2013 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef SPMICKERNEL_HPP
#define SPMICKERNEL_HPP

#ifdef USEMIC
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
#include <iostream>
#include <sstream>
#include <cstring>
#include "parallel/tools/PartitioningTool.hpp"
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

#include "base/grid/GridStorage.hpp"
#include "parallel/datadriven/basis/common/mic/MICKernelBase.hpp"

#include "parallel/datadriven/basis/common/mic/SPMICKernelImpl.hpp"

#ifdef X86_MIC_SYMMETRIC
#include "parallel/datadriven/basis/common/SPX86SimdKernelBase.hpp"
#include "base/exception/operation_exception.hpp"
#endif

namespace sg {
  namespace parallel {
    template<typename KernelImplementation>
    class SPMICKernel {
      public:
        SPMICKernel() {
          // do initialization of device here. Otherwise, it would be done
          // in the first offload in mult, which is already timed => measurements are biased
          // alternative: set environment variable to OFFLOAD_INIT=on_start
#ifdef __INTEL_OFFLOAD
#pragma offload_transfer target(mic)
#else
#ifdef X86_MIC_SYMMETRIC

          if ((KernelImplementation::getChunkDataPoints() % SPX86SimdKernelBase::getChunkDataPoints()) != 0) {
            std::stringstream s;
            s << "For MIC symmetric execution: the chunksize of the MIC Kernel (" << KernelImplementation::getChunkDataPoints() << ") needs to be a multiple of the chunksize of the CPU Kernel (" << SPX86SimdKernelBase::getChunkDataPoints() << ")!" << std::endl;
            std::cerr << s.str() << std::endl;
            throw sg::base::operation_exception(s.str().c_str());
          }

#endif
#endif
          std::cout << "Chunksize of MIC Kernel is " << KernelImplementation::getChunkDataPoints() << std::endl;

          if (micsp::multicard_multtrans_fast) {
            micsp::tempgrid = new float*[micsp::number_mic_devices];

            for (int i = 0; i < micsp::number_mic_devices; i++) {
              micsp::tempgrid[i] = NULL;
            }
          }
        }

        ~SPMICKernel() {
          micsp::deleteGrid();
          micsp::deleteData();

          if (micsp::multicard_multtrans_fast) {
            for (int i = 0; i < micsp::number_mic_devices; i++) {
              if (micsp::tempgrid[i] != NULL) {
                delete[] micsp::tempgrid[i];
              }
            }

            delete[] micsp::tempgrid;
          }
        }

        static const KernelType kernelType = KernelImplementation::kernelType;
        static inline void mult(
          sg::base::DataMatrixSP* level,
          sg::base::DataMatrixSP* index,
          sg::base::DataMatrixSP* mask,
          sg::base::DataMatrixSP* offset,
          sg::base::DataMatrixSP* dataset,
          sg::base::DataVectorSP& alpha,
          sg::base::DataVectorSP& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {
          size_t result_size = result.getSize();
          size_t dims = dataset->getNrows();
          size_t storageSize = alpha.getSize();
          #pragma omp single
          {
            if (micsp::ptrData == NULL) {
              micsp::uploadData(dataset);
            }

            if (micsp::ptrLevel == NULL) {
              micsp::uploadGrid(level, index, mask, offset);
            }
          }
#ifdef __INTEL_OFFLOAD
          //  multi device offload (done in parallel via OpenMP on host)
          // WARNING only works efficiently if there are more threads on the host than coprocessors in the system atm
          size_t start;
          size_t end;
          PartitioningTool::getOpenMPPartitionSegment(0, micsp::number_mic_devices, &start, &end, 1);

          for (size_t d = start; d < end; d++) {
            size_t start_index_data_micdev;
            size_t end_index_data_micdev;
            PartitioningTool::getPartitionSegment(start_index_data, end_index_data, micsp::number_mic_devices, d, &start_index_data_micdev, &end_index_data_micdev, KernelImplementation::getChunkDataPoints());
            size_t mychunk = end_index_data_micdev - start_index_data_micdev;
            // transfer only the data that is needed for this coprocessor
            micsp::transferInputMult(start_index_grid, end_index_grid - start_index_grid, alpha.getPointer(),
                                     start_index_data_micdev, end_index_data_micdev - start_index_data_micdev, result.getPointer(), d);
#pragma offload target(mic:d) \
  nocopy(micsp::ptrLevel) \
  nocopy(micsp::ptrIndex) \
  nocopy(micsp::ptrMask) \
  nocopy(micsp::ptrOffset) \
  nocopy(micsp::ptrData) \
  nocopy(micsp::ptrAlphaMic) \
  nocopy(micsp::ptrDataMic) \
  in(result_size,dims, start_index_grid, end_index_grid, start_index_data_micdev, end_index_data_micdev)
            {
              // set number of thread via env variables (in mpiexec.hydra call via -genv)
              // -genv MIC_ENV_PREFIX MIC -genv MIC_OMP_NUM_THREADS 240 -genv MIC_KMP_AFFINITY granularity=fine,proclist=[1-240],explicit -genv MIC_KMP_BLOCKTIME 0 -genv MIC_KMP_LIBRARY turnaround
              #pragma omp parallel // this parallel section is run on mic, using the hardwarethreads available there
              {
                size_t start_index_data_mic_thread;
                size_t end_index_data_mic_thread;
                PartitioningTool::getOpenMPPartitionSegment(start_index_data_micdev, end_index_data_micdev,
                    &start_index_data_mic_thread, &end_index_data_mic_thread, KernelImplementation::getChunkDataPoints());
                KernelImplementation::multImpl(micsp::ptrLevel, micsp::ptrIndex, micsp::ptrMask, micsp::ptrOffset, micsp::ptrData, micsp::ptrAlphaMic, micsp::ptrDataMic,
                                               result_size, dims, start_index_grid, end_index_grid, start_index_data_mic_thread, end_index_data_mic_thread);
              }
            }
            micsp::transferResultMult(start_index_data_micdev, mychunk, d, result.getPointer()); // workaround intel compiler bug
          }

#else
          size_t start_index_data_mic_thread;
          size_t end_index_data_mic_thread;
          PartitioningTool::getOpenMPPartitionSegment(start_index_data, end_index_data,
              &start_index_data_mic_thread, &end_index_data_mic_thread, KernelImplementation::getChunkDataPoints());
          KernelImplementation::multImpl(micsp::ptrLevel, micsp::ptrIndex, micsp::ptrMask, micsp::ptrOffset, micsp::ptrData, alpha.getPointer(), result.getPointer(),
                                         result_size, dims, start_index_grid, end_index_grid, start_index_data_mic_thread, end_index_data_mic_thread);
#endif
        }


        static inline void multTranspose(
          sg::base::DataMatrixSP* level,
          sg::base::DataMatrixSP* index,
          sg::base::DataMatrixSP* mask,
          sg::base::DataMatrixSP* offset,
          sg::base::DataMatrixSP* dataset,
          sg::base::DataVectorSP& source,
          sg::base::DataVectorSP& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {
          size_t source_size = source.getSize();
          size_t storageSize = result.getSize();
          size_t dims = dataset->getNrows();

          #pragma omp single
          {
            if (micsp::ptrData == NULL) {
              micsp::uploadData(dataset);
            }

            if (micsp::ptrLevel == NULL) {
              micsp::uploadGrid(level, index, mask, offset);

              if (micsp::multicard_multtrans_fast) {
                for (int i = 0; i < micsp::number_mic_devices; i++) {
                  micsp::tempgrid[i] = new float[result.getSize()];
                }
              }
            }

          }
#ifdef __INTEL_OFFLOAD
          size_t start;
          size_t end;
          PartitioningTool::getOpenMPPartitionSegment(0, micsp::number_mic_devices, &start, &end, 1);

          for (size_t d = start; d < end; d++) { //  multi device offload (done in parallel via OpenMP on host)
            size_t start_index_grid_micdev;
            size_t end_index_grid_micdev;
            size_t start_index_data_micdev;
            size_t end_index_data_micdev;

            if (micsp::multicard_multtrans_fast) {
              start_index_grid_micdev = start_index_grid;
              end_index_grid_micdev = end_index_grid;
              PartitioningTool::getPartitionSegment(start_index_data, end_index_data, micsp::number_mic_devices, d,
                                                    &start_index_data_micdev, &end_index_data_micdev, KernelImplementation::getChunkDataPoints());
            } else {
              start_index_data_micdev = start_index_data;
              end_index_data_micdev = end_index_data;

              PartitioningTool::getPartitionSegment(start_index_grid, end_index_grid, micsp::number_mic_devices, d,
                                                    &start_index_grid_micdev, &end_index_grid_micdev, 1);
            }

            // transfer only the data that is needed for this coprocessor
            micsp::transferInputMultTrans(start_index_data_micdev, end_index_data_micdev - start_index_data_micdev, source.getPointer(),
                                          start_index_grid_micdev, end_index_grid_micdev - start_index_grid_micdev, result.getPointer(), d);
#pragma offload target(mic:d) \
  nocopy(micsp::ptrLevel) \
  nocopy(micsp::ptrIndex) \
  nocopy(micsp::ptrMask) \
  nocopy(micsp::ptrOffset) \
  nocopy(micsp::ptrData) \
  nocopy(micsp::ptrDataMic) \
  nocopy(micsp::ptrAlphaMic) \
  in(source_size,dims, start_index_grid_micdev, end_index_grid_micdev, start_index_data_micdev, end_index_data_micdev, d)
            {
              // set number of thread via env variables (in mpiexec.hydra call via -genv)
              // -genv MIC_ENV_PREFIX MIC -genv MIC_OMP_NUM_THREADS 240 -genv MIC_KMP_AFFINITY granularity=fine,proclist=[1-240],explicit -genv MIC_KMP_BLOCKTIME 0 -genv MIC_KMP_LIBRARY turnaround
              #pragma omp parallel // this parallel section is run on mic, using the hardwarethreads available there
              {
                size_t start_index_grid_mic_thread;
                size_t end_index_grid_mic_thread;
                PartitioningTool::getOpenMPPartitionSegment(start_index_grid_micdev, end_index_grid_micdev,
                    &start_index_grid_mic_thread, &end_index_grid_mic_thread, 1);
                KernelImplementation::multTransposeImpl(micsp::ptrLevel, micsp::ptrIndex, micsp::ptrMask, micsp::ptrOffset, micsp::ptrData, micsp::ptrDataMic, micsp::ptrAlphaMic,
                                                        source_size, dims, start_index_grid_mic_thread, end_index_grid_mic_thread, start_index_data_micdev, end_index_data_micdev);
              }
            }

            if (micsp::multicard_multtrans_fast) {
              if (d > 0) {
                micsp::transferResultMultTrans(start_index_grid_micdev, end_index_grid_micdev - start_index_grid_micdev, d, micsp::tempgrid[d]);
              } else {
                micsp::transferResultMultTrans(start_index_grid_micdev, end_index_grid_micdev - start_index_grid_micdev, d, result.getPointer());
              }
            } else {
              micsp::transferResultMultTrans(start_index_grid_micdev, end_index_grid_micdev - start_index_grid_micdev, d, result.getPointer());
            }

          }

          #pragma omp barrier

          if (micsp::multicard_multtrans_fast) {
            size_t start_index_grid_copy;
            size_t end_index_grid_copy;
            PartitioningTool::getOpenMPPartitionSegment(start_index_grid, end_index_grid, &start_index_grid_copy, &end_index_grid_copy, 1);

            for (size_t i = start_index_grid_copy; i < end_index_grid_copy; i++) {
              for (size_t d = 1; d < micsp::number_mic_devices; d++) {
                result[i] += micsp::tempgrid[d][i];
              }
            }
          }

#else
          size_t start_index_grid_mic_thread;
          size_t end_index_grid_mic_thread;
          PartitioningTool::getOpenMPPartitionSegment(start_index_grid, end_index_grid,
              &start_index_grid_mic_thread, &end_index_grid_mic_thread, 1);
          KernelImplementation::multTransposeImpl(micsp::ptrLevel, micsp::ptrIndex, micsp::ptrMask, micsp::ptrOffset, micsp::ptrData, source.getPointer(), result.getPointer(),
                                                  source_size, dims, start_index_grid_mic_thread, end_index_grid_mic_thread, start_index_data, end_index_data);

#endif
        }
        inline void resetKernel() {
          micsp::deleteGrid();

          if (micsp::multicard_multtrans_fast) {
            for (int i = 0; i < micsp::number_mic_devices; i++) {
              if (micsp::tempgrid[i] != NULL) {
                delete[] micsp::tempgrid[i];
                micsp::tempgrid[i] = NULL;
              }
            }
          }
        }
    };

  }
}
#endif // USEMIC
#endif // SPMICKERNEL_HPP

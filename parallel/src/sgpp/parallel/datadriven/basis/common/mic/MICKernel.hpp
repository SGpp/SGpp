// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MICKERNEL_HPP
#define MICKERNEL_HPP

#ifdef USEMIC
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
#include <iostream>
#include <sstream>
#include <cstring>
#include <sgpp/parallel/tools/PartitioningTool.hpp>
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/parallel/datadriven/basis/common/mic/MICKernelBase.hpp>

#include <sgpp/parallel/datadriven/basis/common/mic/MICKernelImpl.hpp>

#ifdef X86_MIC_SYMMETRIC
#include <sgpp/parallel/datadriven/basis/common/X86SimdKernelBase.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#endif

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace parallel {
template<typename KernelImplementation>
class MICKernel {
 public:
  MICKernel() {
    // do initialization of device here. Otherwise, it would be done
    // in the first offload in mult, which is already timed => measurements are biased
    // alternative: set environment variable to OFFLOAD_INIT=on_start
#ifdef __INTEL_OFFLOAD
#pragma offload_transfer target(mic)
#else
#ifdef X86_MIC_SYMMETRIC

    if ((KernelImplementation::getChunkDataPoints() %
         X86SimdKernelBase::getChunkDataPoints()) != 0) {
      std::stringstream s;
      s << "For MIC symmetric execution: the chunksize of the MIC Kernel (" <<
        KernelImplementation::getChunkDataPoints() <<
        ") needs to be a multiple of the chunksize of the CPU Kernel (" <<
        X86SimdKernelBase::getChunkDataPoints() << ")!";
      throw SGPP::base::operation_exception(s.str().c_str());
    }

#endif
#endif
    std::cout << "Chunksize of MIC Kernel is " <<
              KernelImplementation::getChunkDataPoints() << std::endl;

    if (mic::multicard_multtrans_fast) {
      mic::tempgrid = new double*[mic::number_mic_devices];

      for (int i = 0; i < mic::number_mic_devices; i++) {
        mic::tempgrid[i] = NULL;
      }
    }
  }


  ~MICKernel() {
    mic::deleteGrid();
    mic::deleteData();

    if (mic::multicard_multtrans_fast) {
      for (int i = 0; i < mic::number_mic_devices; i++) {
        if (mic::tempgrid[i] != NULL) {
          delete[] mic::tempgrid[i];
        }
      }

      delete[] mic::tempgrid;
    }
  }

  static const KernelType kernelType = KernelImplementation::kernelType;
  static inline void mult(
    SGPP::base::DataMatrix* level,
    SGPP::base::DataMatrix* index,
    SGPP::base::DataMatrix* mask,
    SGPP::base::DataMatrix* offset,
    SGPP::base::DataMatrix* dataset,
    SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result,
    const size_t start_index_grid,
    const size_t end_index_grid,
    const size_t start_index_data,
    const size_t end_index_data) {
    size_t result_size = result.getSize();
    size_t dims = dataset->getNrows();
    size_t storageSize = alpha.getSize();
    #pragma omp single
    {
      if (mic::ptrData == NULL) {
        mic::uploadData(dataset);
      }

      if (mic::ptrLevel == NULL) {
        mic::uploadGrid(level, index, mask, offset);
      }
    }
#ifdef __INTEL_OFFLOAD
    //  multi device offload (done in parallel via OpenMP on host)
    // WARNING only works efficiently if there are more threads on the host than coprocessors in the system atm
    size_t start;
    size_t end;
    PartitioningTool::getOpenMPPartitionSegment(0, mic::number_mic_devices, &start,
        &end, 1);

    for (size_t d = start; d < end; d++) {
      size_t start_index_data_micdev;
      size_t end_index_data_micdev;
      PartitioningTool::getPartitionSegment(start_index_data, end_index_data,
                                            mic::number_mic_devices, d, &start_index_data_micdev, &end_index_data_micdev,
                                            KernelImplementation::getChunkDataPoints());
      size_t mychunk = end_index_data_micdev - start_index_data_micdev;
      // transfer only the data that is needed for this coprocessor
      mic::transferInputMult(start_index_grid, end_index_grid - start_index_grid,
                             alpha.getPointer(),
                             start_index_data_micdev, end_index_data_micdev - start_index_data_micdev,
                             result.getPointer(), d);
#pragma offload target(mic:d) \
  nocopy(mic::ptrLevel) \
  nocopy(mic::ptrIndex) \
  nocopy(mic::ptrMask) \
  nocopy(mic::ptrOffset) \
  nocopy(mic::ptrData) \
  nocopy(mic::ptrAlphaMic) \
  nocopy(mic::ptrDataMic) \
  in(result_size,dims, start_index_grid, end_index_grid, start_index_data_micdev, end_index_data_micdev)
      {
        // set number of thread via env variables (in mpiexec.hydra call via -genv)
        // -genv MIC_ENV_PREFIX MIC -genv MIC_OMP_NUM_THREADS 240 -genv MIC_KMP_AFFINITY granularity=fine,proclist=[1-240],explicit -genv MIC_KMP_BLOCKTIME 0 -genv MIC_KMP_LIBRARY turnaround
        #pragma omp parallel // this parallel section is run on mic, using the hardwarethreads available there
        {
          size_t start_index_data_mic_thread;
          size_t end_index_data_mic_thread;
          PartitioningTool::getOpenMPPartitionSegment(start_index_data_micdev,
              end_index_data_micdev,
              &start_index_data_mic_thread, &end_index_data_mic_thread,
              KernelImplementation::getChunkDataPoints());
          KernelImplementation::multImpl(mic::ptrLevel, mic::ptrIndex, mic::ptrMask,
                                         mic::ptrOffset, mic::ptrData, mic::ptrAlphaMic, mic::ptrDataMic,
                                         result_size, dims, start_index_grid, end_index_grid,
                                         start_index_data_mic_thread, end_index_data_mic_thread);
        }
      }
      mic::transferResultMult(start_index_data_micdev, mychunk, d,
                              result.getPointer()); // workaround intel compiler bug
    }

#else
    size_t start_index_data_mic_thread;
    size_t end_index_data_mic_thread;
    PartitioningTool::getOpenMPPartitionSegment(start_index_data, end_index_data,
        &start_index_data_mic_thread, &end_index_data_mic_thread,
        KernelImplementation::getChunkDataPoints());
    KernelImplementation::multImpl(mic::ptrLevel, mic::ptrIndex, mic::ptrMask,
                                   mic::ptrOffset, mic::ptrData, alpha.getPointer(), result.getPointer(),
                                   result_size, dims, start_index_grid, end_index_grid,
                                   start_index_data_mic_thread, end_index_data_mic_thread);
#endif
  }


  static inline void multTranspose(
    SGPP::base::DataMatrix* level,
    SGPP::base::DataMatrix* index,
    SGPP::base::DataMatrix* mask,
    SGPP::base::DataMatrix* offset,
    SGPP::base::DataMatrix* dataset,
    SGPP::base::DataVector& source,
    SGPP::base::DataVector& result,
    const size_t start_index_grid,
    const size_t end_index_grid,
    const size_t start_index_data,
    const size_t end_index_data) {
    size_t source_size = source.getSize();
    size_t storageSize = result.getSize();
    size_t dims = dataset->getNrows();

    #pragma omp single
    {
      if (mic::ptrData == NULL) {
        mic::uploadData(dataset);
      }

      if (mic::ptrLevel == NULL) {
        mic::uploadGrid(level, index, mask, offset);

        if (mic::multicard_multtrans_fast) {
          for (int i = 0; i < mic::number_mic_devices; i++) {
            mic::tempgrid[i] = new double[result.getSize()];
          }
        }
      }

    }
#ifdef __INTEL_OFFLOAD
    size_t start;
    size_t end;
    PartitioningTool::getOpenMPPartitionSegment(0, mic::number_mic_devices, &start,
        &end, 1);

    for (size_t d = start; d < end;
         d++) { //  multi device offload (done in parallel via OpenMP on host)
      size_t start_index_grid_micdev;
      size_t end_index_grid_micdev;
      size_t start_index_data_micdev;
      size_t end_index_data_micdev;

      if (mic::multicard_multtrans_fast) {
        start_index_grid_micdev = start_index_grid;
        end_index_grid_micdev = end_index_grid;
        PartitioningTool::getPartitionSegment(start_index_data, end_index_data,
                                              mic::number_mic_devices, d,
                                              &start_index_data_micdev, &end_index_data_micdev,
                                              KernelImplementation::getChunkDataPoints());
      } else {
        start_index_data_micdev = start_index_data;
        end_index_data_micdev = end_index_data;

        PartitioningTool::getPartitionSegment(start_index_grid, end_index_grid,
                                              mic::number_mic_devices, d,
                                              &start_index_grid_micdev, &end_index_grid_micdev, 1);
      }

      // transfer only the data that is needed for this coprocessor
      mic::transferInputMultTrans(start_index_data_micdev,
                                  end_index_data_micdev - start_index_data_micdev, source.getPointer(),
                                  start_index_grid_micdev, end_index_grid_micdev - start_index_grid_micdev,
                                  result.getPointer(), d);
#pragma offload target(mic:d) \
  nocopy(mic::ptrLevel) \
  nocopy(mic::ptrIndex) \
  nocopy(mic::ptrMask) \
  nocopy(mic::ptrOffset) \
  nocopy(mic::ptrData) \
  nocopy(mic::ptrDataMic) \
  nocopy(mic::ptrAlphaMic) \
  in(source_size,dims, start_index_grid_micdev, end_index_grid_micdev, start_index_data_micdev, end_index_data_micdev, d)
      {
        // set number of thread via env variables (in mpiexec.hydra call via -genv)
        // -genv MIC_ENV_PREFIX MIC -genv MIC_OMP_NUM_THREADS 240 -genv MIC_KMP_AFFINITY granularity=fine,proclist=[1-240],explicit -genv MIC_KMP_BLOCKTIME 0 -genv MIC_KMP_LIBRARY turnaround
        #pragma omp parallel // this parallel section is run on mic, using the hardwarethreads available there
        {
          size_t start_index_grid_mic_thread;
          size_t end_index_grid_mic_thread;
          PartitioningTool::getOpenMPPartitionSegment(start_index_grid_micdev,
              end_index_grid_micdev,
              &start_index_grid_mic_thread, &end_index_grid_mic_thread, 1);
          KernelImplementation::multTransposeImpl(mic::ptrLevel, mic::ptrIndex,
                                                  mic::ptrMask, mic::ptrOffset, mic::ptrData, mic::ptrDataMic, mic::ptrAlphaMic,
                                                  source_size, dims, start_index_grid_mic_thread, end_index_grid_mic_thread,
                                                  start_index_data_micdev, end_index_data_micdev);
        }
      }

      if (mic::multicard_multtrans_fast) {
        if (d > 0) {
          mic::transferResultMultTrans(start_index_grid_micdev,
                                       end_index_grid_micdev - start_index_grid_micdev, d, mic::tempgrid[d]);
        } else {
          mic::transferResultMultTrans(start_index_grid_micdev,
                                       end_index_grid_micdev - start_index_grid_micdev, d, result.getPointer());
        }
      } else {
        mic::transferResultMultTrans(start_index_grid_micdev,
                                     end_index_grid_micdev - start_index_grid_micdev, d, result.getPointer());
      }

    }

    #pragma omp barrier

    if (mic::multicard_multtrans_fast) {
      size_t start_index_grid_copy;
      size_t end_index_grid_copy;
      PartitioningTool::getOpenMPPartitionSegment(start_index_grid, end_index_grid,
          &start_index_grid_copy, &end_index_grid_copy, 1);

      for (size_t i = start_index_grid_copy; i < end_index_grid_copy; i++) {
        for (size_t d = 1; d < mic::number_mic_devices; d++) {
          result[i] += mic::tempgrid[d][i];
        }
      }
    }

#else
    size_t start_index_grid_mic_thread;
    size_t end_index_grid_mic_thread;
    PartitioningTool::getOpenMPPartitionSegment(start_index_grid, end_index_grid,
        &start_index_grid_mic_thread, &end_index_grid_mic_thread, 1);
    KernelImplementation::multTransposeImpl(mic::ptrLevel, mic::ptrIndex,
                                            mic::ptrMask, mic::ptrOffset, mic::ptrData, source.getPointer(),
                                            result.getPointer(),
                                            source_size, dims, start_index_grid_mic_thread, end_index_grid_mic_thread,
                                            start_index_data, end_index_data);

#endif
  }
  inline void resetKernel() {
    mic::deleteGrid();

    if (mic::multicard_multtrans_fast) {
      for (int i = 0; i < mic::number_mic_devices; i++) {
        if (mic::tempgrid[i] != NULL) {
          delete[] mic::tempgrid[i];
          mic::tempgrid[i] = NULL;
        }
      }
    }
  }
};

}
}
#endif // USEMIC
#endif // MICKERNEL_HPP
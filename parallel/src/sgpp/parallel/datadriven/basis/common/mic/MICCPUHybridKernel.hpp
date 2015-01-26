/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef MICCPUHYBRIDKERNEL_HPP
#define MICCPUHYBRIDKERNEL_HPP

// MIC is required
#ifdef USEMIC

// only makes sense in offload mode
#ifdef __INTEL_OFFLOAD

// OpenMP is required for hybrid execution
#ifdef _OPENMP
#include <omp.h>

#include "parallel/tools/DynamicTwoPartitionAutoTuning.hpp"
#include "base/tools/SGppStopwatch.hpp"
#include "parallel/datadriven/basis/common/mic/MICKernelImpl.hpp"
#include "base/exception/operation_exception.hpp"

namespace sg {
  namespace parallel {

    template<typename CPUImplementation, typename MICImplementation>
    class MICCPUHybridKernel {
      public:
        MICCPUHybridKernel() {
          if ( (MICImplementation::getChunkDataPoints() % CPUImplementation::getChunkDataPoints()) != 0) {
            throw sg::base::operation_exception("For MIC Hybrid Kernel: the chunksize of the MIC Kernel needs to be a multiple of the chunksize of the CPU Kernel!");
          }

          std::cout << "Chunksize of MIC Kernel is " << MICImplementation::getChunkDataPoints() << std::endl;
          tuningMult = new sg::parallel::DynamicTwoPartitionAutoTuning(1, MICImplementation::getChunkDataPoints(), 10);
          tuningMultTrans = new sg::parallel::DynamicTwoPartitionAutoTuning(1, 1, 10);;
          int num_threads = 1;
          #pragma omp parallel
          {
            num_threads = omp_get_num_threads();
          }

          if (num_threads == 1) {
            std::cout << "MIC Hybrid Kernel needs to be executed with at least two threads." << std::endl;
            throw sg::base::operation_exception("MIC Hybrid Kernel needs to be executed with at least two threads.");
          }

          cpu_times = new double[num_threads];
        }
        ~MICCPUHybridKernel() {
          delete[] cpu_times;
          mic::deleteGrid();
          mic::deleteData();
          delete tuningMult;
          delete tuningMultTrans;
        }

        static const KernelType kernelType = CPUImplementation::kernelType; // should match MICImplementation::kernelType
        inline void mult(
          sg::base::DataMatrix* level,
          sg::base::DataMatrix* index,
          sg::base::DataMatrix* mask,
          sg::base::DataMatrix* offset,
          sg::base::DataMatrix* dataset,
          sg::base::DataVector& alpha,
          sg::base::DataVector& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {
          int tid = omp_get_thread_num();
          int num_threads = omp_get_num_threads();
          size_t result_size = result.getSize();
          size_t dims = dataset->getNrows();
          size_t data_range = end_index_data - start_index_data;

          // initialize loadbalancing
          if (tuningMult->getProblemSize() != data_range) {
            #pragma omp barrier
            #pragma omp master
            {
              tuningMult->setProblemSize(data_range);
            }
            #pragma omp barrier
          }

          double mic_time = 0.0;

          cpu_times[tid] = 0.0;

          // split result into MIC and CPU partition
          size_t mic_size = data_range - tuningMult->getPartition1Size();
          size_t end_index_data_mic = start_index_data + mic_size;

          if (tid == 0) {
            if (mic::ptrData == NULL) {
              mic::uploadData(dataset);
            }

            if (mic::ptrLevel == NULL) {
              mic::uploadGrid(level, index, mask, offset);
            }


            double mic_start = omp_get_wtime();
            size_t* start_indexes_data_micdev = new size_t[mic::number_mic_devices];
            size_t* end_indexes_data_micdev = new size_t[mic::number_mic_devices];
            int* finished = new int[mic::number_mic_devices];

            for (size_t d = 0; d < mic::number_mic_devices; d++) { // async multi device offload
              PartitioningTool::getPartitionSegment(start_index_data, end_index_data_mic, mic::number_mic_devices, d, &(start_indexes_data_micdev[d]), &(end_indexes_data_micdev[d]), MICImplementation::getChunkDataPoints());
              size_t start_index_data_micdev = start_indexes_data_micdev[d];
              size_t end_index_data_micdev = end_indexes_data_micdev[d];
              // transfer only the data that is needed for this coprocessor
              mic::transferInputMult(start_index_grid, end_index_grid - start_index_grid, alpha.getPointer(),
                                     start_index_data_micdev, end_index_data_micdev - start_index_data_micdev, result.getPointer(), d);
#pragma offload target(mic:d) \
    nocopy(mic::ptrLevel) \
    nocopy(mic::ptrIndex) \
    nocopy(mic::ptrMask) \
    nocopy(mic::ptrOffset) \
    nocopy(mic::ptrData) \
    nocopy(mic::ptrAlphaMic) \
    nocopy(mic::ptrDataMic) \
    in(result_size,dims, start_index_grid, end_index_grid, start_index_data_micdev, end_index_data_micdev) \
    signal(&(finished[d]))
              {
                // set number of thread via env variables (in mpiexec.hydra call via -genv)
                // -genv MIC_ENV_PREFIX MIC -genv MIC_OMP_NUM_THREADS 240 -genv MIC_KMP_AFFINITY granularity=fine,proclist=[1-240],explicit -genv MIC_KMP_BLOCKTIME 0 -genv MIC_KMP_LIBRARY turnaround
                #pragma omp parallel // this parallel section is run on mic, using the hardwarethreads available there
                {
                  size_t start_index_data_mic_thread;
                  size_t end_index_data_mic_thread;
                  PartitioningTool::getOpenMPPartitionSegment(start_index_data_micdev, end_index_data_micdev,
                      &start_index_data_mic_thread, &end_index_data_mic_thread, MICImplementation::getChunkDataPoints());
                  MICImplementation::multImpl(mic::ptrLevel, mic::ptrIndex, mic::ptrMask, mic::ptrOffset, mic::ptrData, mic::ptrAlphaMic, mic::ptrDataMic,
                                              result_size, dims, start_index_grid, end_index_grid, start_index_data_mic_thread, end_index_data_mic_thread);
                }
              }
            }

            for (size_t d = 0; d < mic::number_mic_devices; d++) {
#pragma offload_wait target(mic:d) wait(&(finished[d]))
              mic::transferResultMult(start_indexes_data_micdev[d], end_indexes_data_micdev[d] - start_indexes_data_micdev[d], d, result.getPointer()); // workaround intel compiler bug
            }

            mic_time = omp_get_wtime() - mic_start;
            delete[] start_indexes_data_micdev;
            delete[] end_indexes_data_micdev;
            delete[] finished;
          } else {
            // distribute work evenly across all threads but thread 0 (it handles OpenCL kernels)
            size_t cpu_size = end_index_data - end_index_data_mic;

            if ((cpu_size % CPUImplementation::getChunkDataPoints()) != 0) {
              throw sg::base::operation_exception("MIC Blocksize must be multiple of CPU Blocksize for Hybrid to function correctly");
            }

            size_t start_index_data_x86simd = end_index_data_mic;
            size_t end_index_data_x86simd = end_index_data;

            size_t start_index_data_thread;
            size_t end_index_data_thread;

            PartitioningTool::getPartitionSegment(start_index_data_x86simd, end_index_data_x86simd, num_threads - 1, tid - 1, &start_index_data_thread, &end_index_data_thread, CPUImplementation::getChunkDataPoints());
            double start = omp_get_wtime();

            CPUImplementation::multImpl(level, index, mask, offset, dataset, alpha, result,
                                        start_index_grid, end_index_grid, start_index_data_thread, end_index_data_thread);
            cpu_times[tid] = omp_get_wtime() - start;
          }

          #pragma omp barrier
          #pragma omp master
          {
            double max_cpu_time = 0.0;

            for (int i = 0; i < num_threads; i++) {
              if (cpu_times[i] > max_cpu_time) {
                max_cpu_time = cpu_times[i];
              }
            }

            tuningMult->setExecutionTimes(max_cpu_time, mic_time);
          }
          #pragma omp barrier
        }

        inline void multTranspose(
          sg::base::DataMatrix* level,
          sg::base::DataMatrix* index,
          sg::base::DataMatrix* mask,
          sg::base::DataMatrix* offset,
          sg::base::DataMatrix* dataset,
          sg::base::DataVector& source,
          sg::base::DataVector& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {
          size_t grid_range = end_index_grid - start_index_grid;
          size_t source_size = source.getSize();
          size_t dims = dataset->getNrows();

          int tid = omp_get_thread_num();
          int num_threads = omp_get_num_threads();

          double mic_time = 0.0;
          cpu_times[tid] = 0.0;

          if (tuningMultTrans->getProblemSize() != grid_range) {
            #pragma omp barrier
            #pragma omp master
            {
              tuningMultTrans->setProblemSize(grid_range);
            }
            #pragma omp barrier
          }

          // split result into MIC and CPU partition
          size_t mic_size = grid_range - tuningMultTrans->getPartition1Size();
          size_t end_index_grid_mic = start_index_grid + mic_size;

          if (tid == 0) {
            if (mic::ptrData == NULL) {
              mic::uploadData(dataset);
            }

            if (mic::ptrLevel == NULL) {
              mic::uploadGrid(level, index, mask, offset);
            }

            double mic_start = omp_get_wtime();
            size_t* start_indexes_grid_micdev = new size_t[mic::number_mic_devices];
            size_t* end_indexes_grid_micdev = new size_t[mic::number_mic_devices];
            double* finished = new double[mic::number_mic_devices];

            for (size_t d = 0; d < mic::number_mic_devices; d++) { // async multi device offload
              PartitioningTool::getPartitionSegment(start_index_grid, end_index_grid_mic, mic::number_mic_devices, d, &(start_indexes_grid_micdev[d]), &(end_indexes_grid_micdev[d]), 1);
              size_t start_index_grid_micdev = start_indexes_grid_micdev[d];
              size_t end_index_grid_micdev = end_indexes_grid_micdev[d];
              // transfer only the data that is needed for this coprocessor
              mic::transferInputMultTrans(start_index_data, end_index_data - start_index_data, source.getPointer(),
                                          start_index_grid_micdev, end_index_grid_micdev - start_index_grid_micdev, result.getPointer(), d);
#pragma offload target(mic:d) \
    nocopy(mic::ptrLevel) \
    nocopy(mic::ptrIndex) \
    nocopy(mic::ptrMask) \
    nocopy(mic::ptrOffset) \
    nocopy(mic::ptrData) \
    nocopy(mic::ptrDataMic) \
    nocopy(mic::ptrAlphaMic) \
    in(source_size,dims, start_index_grid_micdev, end_index_grid_micdev, start_index_data, end_index_data, d) \
    signal(&(finished[d]))
              {
                // set number of thread via env variables (in mpiexec.hydra call via -genv)
                // -genv MIC_ENV_PREFIX MIC -genv MIC_OMP_NUM_THREADS 240 -genv MIC_KMP_AFFINITY granularity=fine,proclist=[1-240],explicit -genv MIC_KMP_BLOCKTIME 0 -genv MIC_KMP_LIBRARY turnaround
                #pragma omp parallel // this parallel section is run on mic, using the hardwarethreads available there
                {
                  size_t start_index_grid_mic_thread;
                  size_t end_index_grid_mic_thread;
                  PartitioningTool::getOpenMPPartitionSegment(start_index_grid_micdev, end_index_grid_micdev,
                      &start_index_grid_mic_thread, &end_index_grid_mic_thread, 1);
                  MICImplementation::multTransposeImpl(mic::ptrLevel, mic::ptrIndex, mic::ptrMask, mic::ptrOffset, mic::ptrData, mic::ptrDataMic, mic::ptrAlphaMic,
                                                       source_size, dims, start_index_grid_mic_thread, end_index_grid_mic_thread, start_index_data, end_index_data);
                }
              }
            }

            for (size_t d = 0; d < mic::number_mic_devices; d++) { //  multi device offload
#pragma offload_wait target(mic:d) wait(&(finished[d]))
              mic::transferResultMultTrans(start_indexes_grid_micdev[d], end_indexes_grid_micdev[d] - start_indexes_grid_micdev[d], d, result.getPointer());
            }

            mic_time = omp_get_wtime() - mic_start;
            delete[] start_indexes_grid_micdev;
            delete[] end_indexes_grid_micdev;
            delete[] finished;
          } else {
            size_t start_grid_thread;
            size_t end_grid_thread;
            // distribute work evenly across all threads but thread 0
            PartitioningTool::getPartitionSegment(end_index_grid_mic, end_index_grid, num_threads - 1, tid - 1, &start_grid_thread, &end_grid_thread, 1);
            double start = omp_get_wtime();
            CPUImplementation::multTransposeImpl(level, index, mask, offset, dataset, source, result,
                                                 start_grid_thread, end_grid_thread, start_index_data, end_index_data);
            cpu_times[tid] = omp_get_wtime() - start;
          }

          #pragma omp barrier
          #pragma omp master
          {
            double max_cpu_time = 0.0;

            for (int i = 0; i < num_threads; i++) {
              if (cpu_times[i] > max_cpu_time) {
                max_cpu_time = cpu_times[i];
              }
            }

            tuningMultTrans->setExecutionTimes(max_cpu_time, mic_time);
          }
          #pragma omp barrier
        }

        inline void resetKernel() {
          mic::deleteGrid();
        }

      private:
        double* cpu_times;
        MICImplementation m_oclkernel;
        sg::base::SGppStopwatch myTimer;

        /// Autotuning object for mult routine
        sg::parallel::TwoPartitionAutoTuning* tuningMult;
        /// Autotuning object for mult trans routine
        sg::parallel::TwoPartitionAutoTuning* tuningMultTrans;
    };
  }
}

#endif // _OPENMP

#endif // USEMIC

#endif // __INTEL_OFFLOAD

#endif // MICCPUHYBRIDKERNEL_HPP

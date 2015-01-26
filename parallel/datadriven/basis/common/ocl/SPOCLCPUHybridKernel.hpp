/* ****************************************************************************
* Copyright (C) 2010-2013 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef SPOCLCPUHYBRIDKERNEL_HPP
#define SPOCLCPUHYBRIDKERNEL_HPP

// OpenMP is required for hybrid execution
#ifdef _OPENMP
#include <omp.h>

#include "base/grid/GridStorage.hpp"
#include "parallel/tools/DynamicTwoPartitionAutoTuning.hpp"
#include "parallel/tools/PartitioningTool.hpp"
#include "parallel/datadriven/basis/common/ocl/SPOCLKernelImpl.hpp"
#include "base/exception/operation_exception.hpp"

namespace sg {
  namespace parallel {

    template<typename CPUImplementation, typename OCLBasisType>
    class SPOCLCPUHybridKernel {
      public:
        SPOCLCPUHybridKernel() {
          tuningMult = new sg::parallel::DynamicTwoPartitionAutoTuning(1, m_oclkernel.getChunkDataPoints(), 10);;
          tuningMultTrans = new sg::parallel::DynamicTwoPartitionAutoTuning(1, m_oclkernel.getChunkGridPoints(), 10);;
          int num_threads = 1;
          #pragma omp parallel
          {
            num_threads = omp_get_num_threads();
          }

          if (num_threads == 1) {
            throw sg::base::operation_exception("OpenCL Hybrid Kernel needs to be executed with at least two threads");
          }

          cpu_times = new double[num_threads];
        }
        ~SPOCLCPUHybridKernel() {
          delete[] cpu_times;
          delete tuningMult;
          delete tuningMultTrans;
        }

        static const KernelType kernelType = CPUImplementation::kernelType;
        inline void mult(
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
          int tid = omp_get_thread_num();
          int num_threads = omp_get_num_threads();
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

          double gpu_time = 0.0;

          cpu_times[tid] = 0.0;

          // split result into GPU and CPU partition
          size_t gpu_size = data_range - tuningMult->getPartition1Size();
          size_t end_index_data_gpu = start_index_data + gpu_size;

          if (tid == 0) {
            double loc_start = omp_get_wtime();
            m_oclkernel.multImpl(level, index, mask, offset, dataset, alpha, result,
                                 start_index_grid, end_index_grid, start_index_data, end_index_data_gpu);
            gpu_time = omp_get_wtime() - loc_start;
          } else {
            // distribute work evenly across all threads but thread 0 (it starts OpenCL kernels)
            size_t cpu_size = end_index_data - end_index_data_gpu;
            size_t blocksize = CPUImplementation::getChunkDataPoints();
            cpu_size = (cpu_size / blocksize ) * blocksize;

            size_t start_index_data_x86simd = end_index_data_gpu;
            size_t end_index_data_x86simd = start_index_data_x86simd + cpu_size;

            size_t start_index_data_thread;
            size_t end_index_data_thread;

            PartitioningTool::getPartitionSegment(start_index_data_x86simd, end_index_data_x86simd, num_threads - 1, tid - 1, &start_index_data_thread, &end_index_data_thread, blocksize);
            double start = omp_get_wtime();

            CPUImplementation::multImpl(level, index, mask, offset, dataset, alpha, result,
                                        start_index_grid, end_index_grid, start_index_data_thread, end_index_data_thread);

            if (tid == num_threads - 1) {
              // handle possible remainder
              OCLBasisType::multDefault(level, index, mask, offset, dataset, alpha, result,
                                        start_index_grid, end_index_grid, end_index_data_x86simd, end_index_data);
            }

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

            tuningMult->setExecutionTimes(max_cpu_time, gpu_time);
          }
          #pragma omp barrier
        }

        inline void multTranspose(
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
          size_t grid_range = end_index_grid - start_index_grid;

          int tid = omp_get_thread_num();
          int num_threads = omp_get_num_threads();

          double gpu_time = 0.0;
          cpu_times[tid] = 0.0;

          if (tuningMultTrans->getProblemSize() != grid_range) {
            #pragma omp barrier
            #pragma omp master
            {
              tuningMultTrans->setProblemSize(grid_range);
            }
            #pragma omp barrier
          }

          // split result into GPU and CPU partition
          size_t gpu_size = grid_range - tuningMultTrans->getPartition1Size();
          size_t end_index_grid_gpu = start_index_grid + gpu_size;

          if (tid == 0) {
            double loc_start = omp_get_wtime();
            m_oclkernel.multTransposeImpl(level, index, mask, offset, dataset, source, result,
                                          start_index_grid, end_index_grid_gpu, start_index_data, end_index_data);
            gpu_time = omp_get_wtime() - loc_start;

          } else {
            size_t start_grid_thread;
            size_t end_grid_thread;
            // distribute work evenly across all threads but thread 0
            PartitioningTool::getPartitionSegment(end_index_grid_gpu, end_index_grid, num_threads - 1, tid - 1, &start_grid_thread, &end_grid_thread, 1);
            double start = omp_get_wtime();

            size_t data_range = end_index_data - start_index_data;
            size_t data_range_vec = (data_range / CPUImplementation::getChunkDataPoints()) * CPUImplementation::getChunkDataPoints();
            CPUImplementation::multTransposeImpl(level, index, mask, offset, dataset, source, result,
                                                 start_grid_thread, end_grid_thread, start_index_data, start_index_data + data_range_vec);
            // handle remainder of dataset
            OCLBasisType::multTransposeDefault(level, index, mask, offset, dataset, source, result,
                                               start_grid_thread, end_grid_thread, start_index_data + data_range_vec, end_index_data);
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

            tuningMultTrans->setExecutionTimes(max_cpu_time, gpu_time);
          }
          #pragma omp barrier
        }

        inline void resetKernel() {
          m_oclkernel.resetKernel();
        }

      private:
        double* cpu_times;
        SPOCLKernelImpl<OCLBasisType> m_oclkernel;
        sg::base::SGppStopwatch myTimer;

        /// Autotuning object for mult routine
        sg::parallel::TwoPartitionAutoTuning* tuningMult;
        /// Autotuning object for mult trans routine
        sg::parallel::TwoPartitionAutoTuning* tuningMultTrans;
    };
  }
}
#endif

#endif // SPOCLCPUHYBRIDKERNEL_HPP

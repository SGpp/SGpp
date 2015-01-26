/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef OCLKERNEL_HPP
#define OCLKERNEL_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>
#include <sgpp/parallel/datadriven/basis/common/ocl/OCLKernelImpl.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sg {
  namespace parallel {

    template<typename OCLBasisType>
    class OCLKernel {
      public:
        static const KernelType kernelType = OCLBasisType::kernelType;
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
          #pragma omp master
          {
            //double time = m_oclkernel.multImpl(level, index, mask, offset, dataset, alpha, result,
            m_oclkernel.multImpl(level, index, mask, offset, dataset, alpha, result,
            start_index_grid, end_index_grid, start_index_data, end_index_data);
          }
        }

        inline void multTranspose(
          sg::base::DataMatrix* level,
          sg::base::DataMatrix* index,
          sg::base::DataMatrix* mask, //unused for this specialization
          sg::base::DataMatrix* offset, //unused for this specialization
          sg::base::DataMatrix* dataset,
          sg::base::DataVector& source,
          sg::base::DataVector& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {

#ifdef _OPENMP
          int tid = omp_get_thread_num();
          int num_threads = omp_get_num_threads();
#else
          int tid = 0;
          int num_threads = 1;
#endif

          size_t range = end_index_grid - start_index_grid;
          size_t numWGs = range / m_oclkernel.getChunkGridPoints();
          size_t end_index_grid_gpu = start_index_grid + numWGs *  m_oclkernel.getChunkGridPoints();

          if (tid == 0) {
            //double time = 0.0;
            //time = m_oclkernel.multTransposeImpl(level, index, mask, offset, dataset, source, result,
            m_oclkernel.multTransposeImpl(level, index, mask, offset, dataset, source, result,
                                          start_index_grid, end_index_grid_gpu, start_index_data, end_index_data);
          }

          size_t start_grid_cpu = end_index_grid_gpu;
          size_t end_grid_cpu = end_index_grid;

          if (num_threads > 1) {
            if (tid == 0) {
              // don't do anything for thread 0 when there is more than one thread
              start_grid_cpu = end_grid_cpu = 0;
            } else {
              // distribute work evenly across all threads but thread 0
              PartitioningTool::getPartitionSegment(end_index_grid_gpu, end_index_grid, num_threads - 1, tid - 1, &start_grid_cpu, &end_grid_cpu, 1);
            }
          }

          OCLBasisType::multTransposeDefault(level, index, mask, offset, dataset, source, result,
                                             start_grid_cpu, end_grid_cpu, start_index_data, end_index_data);
        }

        void resetKernel() {
          m_oclkernel.resetKernel();
        }

        OCLKernelImpl<OCLBasisType> m_oclkernel;
    };
  }
}
#endif // OCLKERNEL_HPP

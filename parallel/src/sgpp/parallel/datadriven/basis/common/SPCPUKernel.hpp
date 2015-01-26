/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef SPCPUKERNEL_HPP
#define SPCPUKERNEL_HPP

#include <sgpp/parallel/datadriven/basis/common/SPX86SimdKernelBase.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>

namespace sg {
  namespace parallel {

    template<typename KernelImplementation>
    class SPCPUKernel {
      public:
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

          size_t start;
          size_t end;
          PartitioningTool::getOpenMPPartitionSegment(start_index_data, end_index_data, &start, &end,
              KernelImplementation::getChunkDataPoints());
          KernelImplementation::multImpl(level, index, mask, offset, dataset, alpha, result,
                                         start_index_grid, end_index_grid, start, end);
        }
        static inline void multTranspose(
          sg::base::DataMatrixSP* level,
          sg::base::DataMatrixSP* index,
          sg::base::DataMatrixSP* mask, //unused for this specialization
          sg::base::DataMatrixSP* offset, //unused for this specialization
          sg::base::DataMatrixSP* dataset,
          sg::base::DataVectorSP& source,
          sg::base::DataVectorSP& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {

          size_t start;
          size_t end;
          PartitioningTool::getOpenMPPartitionSegment(start_index_grid, end_index_grid, &start, &end, 1);
          KernelImplementation::multTransposeImpl(level, index, mask, offset, dataset, source, result,
                                                  start, end, start_index_data, end_index_data);
        }
        static inline void resetKernel() {}
    };
  }
}

#endif // SPCPUKERNEL_HPP

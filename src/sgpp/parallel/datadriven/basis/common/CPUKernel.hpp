/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef CPUKERNEL_HPP
#define CPUKERNEL_HPP

#include "base/grid/GridStorage.hpp"
#include "parallel/tools/PartitioningTool.hpp"

namespace sg {
  namespace parallel {

    template<typename KernelImplementation>
    class CPUKernel {
      public:
        static const KernelType kernelType = KernelImplementation::kernelType;
        static inline void mult(
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

          size_t start;
          size_t end;
          PartitioningTool::getOpenMPPartitionSegment(start_index_data, end_index_data, &start, &end,
              KernelImplementation::getChunkDataPoints());
          KernelImplementation::multImpl(level, index, mask, offset, dataset, alpha, result,
                                         start_index_grid, end_index_grid, start, end);
        }
        static inline void multTranspose(
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


#endif // CPUKERNEL_HPP

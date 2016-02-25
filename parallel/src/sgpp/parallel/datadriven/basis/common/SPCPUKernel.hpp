// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SPCPUKERNEL_HPP
#define SPCPUKERNEL_HPP

#include <sgpp/parallel/datadriven/basis/common/SPX86SimdKernelBase.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace parallel {

template <typename KernelImplementation>
class SPCPUKernel {
 public:
  static const KernelType kernelType = KernelImplementation::kernelType;
  static inline void mult(SGPP::base::DataMatrixSP* level, SGPP::base::DataMatrixSP* index,
                          SGPP::base::DataMatrixSP* mask, SGPP::base::DataMatrixSP* offset,
                          SGPP::base::DataMatrixSP* dataset, SGPP::base::DataVectorSP& alpha,
                          SGPP::base::DataVectorSP& result, const size_t start_index_grid,
                          const size_t end_index_grid, const size_t start_index_data,
                          const size_t end_index_data) {
    size_t start;
    size_t end;
    PartitioningTool::getOpenMPPartitionSegment(start_index_data, end_index_data, &start, &end,
                                                KernelImplementation::getChunkDataPoints());
    KernelImplementation::multImpl(level, index, mask, offset, dataset, alpha, result,
                                   start_index_grid, end_index_grid, start, end);
  }
  static inline void multTranspose(
      SGPP::base::DataMatrixSP* level, SGPP::base::DataMatrixSP* index,
      SGPP::base::DataMatrixSP* mask,    // unused for this specialization
      SGPP::base::DataMatrixSP* offset,  // unused for this specialization
      SGPP::base::DataMatrixSP* dataset, SGPP::base::DataVectorSP& source,
      SGPP::base::DataVectorSP& result, const size_t start_index_grid, const size_t end_index_grid,
      const size_t start_index_data, const size_t end_index_data) {
    size_t start;
    size_t end;
    PartitioningTool::getOpenMPPartitionSegment(start_index_grid, end_index_grid, &start, &end, 1);
    KernelImplementation::multTransposeImpl(level, index, mask, offset, dataset, source, result,
                                            start, end, start_index_data, end_index_data);
  }
  static inline void resetKernel() {}
};
}  // namespace parallel
}  // namespace SGPP

#endif  // SPCPUKERNEL_HPP

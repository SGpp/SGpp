// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SPCPUKERNEL_HPP
#define SPCPUKERNEL_HPP

#include <sgpp/parallel/datadriven/basis/common/SPX86SimdKernelBase.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

template <typename KernelImplementation>
class SPCPUKernel {
 public:
  static const KernelType kernelType = KernelImplementation::kernelType;
  static inline void mult(sgpp::base::DataMatrixSP* level, sgpp::base::DataMatrixSP* index,
                          sgpp::base::DataMatrixSP* mask, sgpp::base::DataMatrixSP* offset,
                          sgpp::base::DataMatrixSP* dataset, sgpp::base::DataVectorSP& alpha,
                          sgpp::base::DataVectorSP& result, const size_t start_index_grid,
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
      sgpp::base::DataMatrixSP* level, sgpp::base::DataMatrixSP* index,
      sgpp::base::DataMatrixSP* mask,    // unused for this specialization
      sgpp::base::DataMatrixSP* offset,  // unused for this specialization
      sgpp::base::DataMatrixSP* dataset, sgpp::base::DataVectorSP& source,
      sgpp::base::DataVectorSP& result, const size_t start_index_grid, const size_t end_index_grid,
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
}  // namespace sgpp

#endif  // SPCPUKERNEL_HPP

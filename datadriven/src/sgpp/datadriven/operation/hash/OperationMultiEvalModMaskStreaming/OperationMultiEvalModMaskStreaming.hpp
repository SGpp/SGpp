// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

#ifndef STREAMING_MODLINEAR_MIC_AVX512_UNROLLING_WIDTH
// #define STREAMING_MODLINEAR_MIC_AVX512_UNROLLING_WIDTH 24
#define STREAMING_MODLINEAR_MIC_AVX512_UNROLLING_WIDTH 96
#endif

namespace sgpp {
namespace datadriven {

class OperationMultiEvalModMaskStreaming : public base::OperationMultipleEval {
 protected:
  sgpp::base::DataMatrix preparedDataset;
  /// Member to store the sparse grid's levels for better vectorization
  std::vector<double> level;
  /// Member to store the sparse grid's indices for better vectorization
  std::vector<double> index;

  std::vector<double> mask;
  std::vector<double> offset;
  /// Timer object to handle time measurements
  sgpp::base::SGppStopwatch myTimer_;

  base::GridStorage* storage;

  double duration;

 public:
  OperationMultiEvalModMaskStreaming(base::Grid& grid,
                                     base::DataMatrix& dataset);

  ~OperationMultiEvalModMaskStreaming() override;

  size_t getChunkGridPoints();

  size_t getChunkDataPoints();

  void mult(sgpp::base::DataVector& alpha,
            sgpp::base::DataVector& result) override;

  void multTranspose(sgpp::base::DataVector& source,
                     sgpp::base::DataVector& result) override;

  void prepare() override;

  double getDuration() override;

 private:
  void getPartitionSegment(size_t start, size_t end, size_t segmentCount,
                           size_t segmentNumber, size_t* segmentStart,
                           size_t* segmentEnd, size_t blockSize);

  size_t padDataset(sgpp::base::DataMatrix& dataset);

  void getOpenMPPartitionSegment(size_t start, size_t end, size_t* segmentStart,
                                 size_t* segmentEnd, size_t blocksize);

  void multImpl(std::vector<double>& level, std::vector<double>& index,
                std::vector<double>& mask, std::vector<double>& offset,
                sgpp::base::DataMatrix* dataset, sgpp::base::DataVector& alpha,
                sgpp::base::DataVector& result, const size_t start_index_grid,
                const size_t end_index_grid, const size_t start_index_data,
                const size_t end_index_data);

  void multTransposeImpl(std::vector<double>& level, std::vector<double>& index,
                         std::vector<double>& mask, std::vector<double>& offset,
                         sgpp::base::DataMatrix* dataset,
                         sgpp::base::DataVector& source,
                         sgpp::base::DataVector& result,
                         const size_t start_index_grid,
                         const size_t end_index_grid,
                         const size_t start_index_data,
                         const size_t end_index_data);

  void recalculateLevelIndexMask();
};
}  // namespace datadriven
}  // namespace sgpp

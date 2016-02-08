// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

#ifndef STREAMING_MODLINEAR_MIC_AVX512_UNROLLING_WIDTH
// #define STREAMING_MODLINEAR_MIC_AVX512_UNROLLING_WIDTH 24
#define STREAMING_MODLINEAR_MIC_AVX512_UNROLLING_WIDTH 96
#endif

namespace SGPP {
namespace datadriven {

class OperationMultiEvalModMaskStreaming: public base::OperationMultipleEval {
 protected:
  SGPP::base::DataMatrix preparedDataset;
  /// Member to store the sparse grid's levels for better vectorization
  std::vector<double> level;
  /// Member to store the sparse grid's indices for better vectorization
  std::vector<double> index;

  std::vector<double> mask;
  std::vector<double> offset;
  /// Timer object to handle time measurements
  SGPP::base::SGppStopwatch myTimer_;

  base::GridStorage* storage;

  float_t duration;

 public:
  OperationMultiEvalModMaskStreaming(base::Grid& grid, base::DataMatrix& dataset);

  ~OperationMultiEvalModMaskStreaming();

  size_t getChunkGridPoints();

  size_t getChunkDataPoints();

  void mult(SGPP::base::DataVector& alpha,
            SGPP::base::DataVector& result) override;

  void multTranspose(SGPP::base::DataVector& source,
                     SGPP::base::DataVector& result) override;

  void prepare() override;

  float_t getDuration() override;

 private:
  void getPartitionSegment(size_t start, size_t end, size_t segmentCount,
                           size_t segmentNumber, size_t* segmentStart,
                           size_t* segmentEnd, size_t blockSize);

  size_t padDataset(SGPP::base::DataMatrix& dataset);

  void getOpenMPPartitionSegment(size_t start, size_t end, size_t* segmentStart,
                                 size_t* segmentEnd,
                                 size_t blocksize);

  void multImpl(std::vector<double>& level, std::vector<double>& index,
                std::vector<double>& mask,
                std::vector<double>& offset,
                SGPP::base::DataMatrix* dataset,
                SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                const size_t start_index_grid,
                const size_t end_index_grid, const size_t start_index_data,
                const size_t end_index_data);

  void multTransposeImpl(std::vector<double>& level, std::vector<double>& index,
                         std::vector<double>& mask,
                         std::vector<double>& offset,
                         SGPP::base::DataMatrix* dataset,
                         SGPP::base::DataVector& source, SGPP::base::DataVector& result,
                         const size_t start_index_grid,
                         const size_t end_index_grid, const size_t start_index_data,
                         const size_t end_index_data);

  void recalculateLevelIndexMask();
};

}  // namespace datadriven
}  // namespace SGPP


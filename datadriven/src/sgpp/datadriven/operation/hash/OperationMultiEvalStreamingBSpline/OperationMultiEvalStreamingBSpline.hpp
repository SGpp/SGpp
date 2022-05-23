#pragma once

#include <omp.h>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/globaldef.hpp>

#ifndef STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH
//#define STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH 24
#define STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH 96
#endif

namespace sgpp {
namespace datadriven {

class OperationMultiEvalStreamingBSpline : public base::OperationMultipleEval {
 protected:
  sgpp::base::DataMatrix preparedDataset;
  /// Member to store the sparse grid's levels for better vectorization
  sgpp::base::DataMatrix* level_ = nullptr;
  /// Member to store the sparse grid's indices for better vectorization
  sgpp::base::DataMatrix* index_ = nullptr;
  /// Timer object to handle time measurements
  sgpp::base::SGppStopwatch myTimer_;

  base::GridStorage* storage;

  double duration;

 public:
  OperationMultiEvalStreamingBSpline(base::Grid& grid, base::DataMatrix& dataset);

  ~OperationMultiEvalStreamingBSpline() override;

  size_t getChunkGridPoints();

  size_t getChunkDataPoints();

  void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) override;

  void multTranspose(sgpp::base::DataVector& source, sgpp::base::DataVector& result) override;

  void prepare() override;

  double getDuration() override;

 private:
  void getPartitionSegment(size_t start, size_t end, size_t segmentCount, size_t segmentNumber,
                           size_t* segmentStart, size_t* segmentEnd, size_t blockSize);

  size_t padDataset(sgpp::base::DataMatrix& dataset);

  void getOpenMPPartitionSegment(size_t start, size_t end, size_t* segmentStart, size_t* segmentEnd,
                                 size_t blocksize);

  void multImpl(sgpp::base::DataMatrix* level, sgpp::base::DataMatrix* index,
                sgpp::base::DataMatrix* dataset, sgpp::base::DataVector& alpha,
                sgpp::base::DataVector& result, const size_t start_index_grid,
                const size_t end_index_grid, const size_t start_index_data,
                const size_t end_index_data);

  void multTransposeImpl(sgpp::base::DataMatrix* level, sgpp::base::DataMatrix* index,
                         sgpp::base::DataMatrix* dataset, sgpp::base::DataVector& source,
                         sgpp::base::DataVector& result, const size_t start_index_grid,
                         const size_t end_index_grid, const size_t start_index_data,
                         const size_t end_index_data);

  void recalculateLevelAndIndex();
};
}  // namespace datadriven
}  // namespace sgpp

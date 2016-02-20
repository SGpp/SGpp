// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultiEvalModMaskStreaming/OperationMultiEvalModMaskStreaming.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace SGPP {
namespace datadriven {

OperationMultiEvalModMaskStreaming::OperationMultiEvalModMaskStreaming(base::Grid& grid,
                                                                       base::DataMatrix& dataset)
    : OperationMultipleEval(grid, dataset),
      preparedDataset(dataset),
      myTimer_(SGPP::base::SGppStopwatch()),
      duration(-1.0) {
  this->storage = &grid.getStorage();
  this->padDataset(this->preparedDataset);
  this->preparedDataset.transpose();

  // create the kernel specific data structures for the current grid
  this->prepare();
}

OperationMultiEvalModMaskStreaming::~OperationMultiEvalModMaskStreaming() {}

void OperationMultiEvalModMaskStreaming::getPartitionSegment(size_t start, size_t end,
                                                             size_t segmentCount,
                                                             size_t segmentNumber,
                                                             size_t* segmentStart,
                                                             size_t* segmentEnd, size_t blockSize) {
  size_t totalSize = end - start;

  // check for valid input
  if (blockSize == 0) {
    throw SGPP::base::operation_exception("blockSize must not be zero!");
  }

  if (totalSize % blockSize != 0) {
    // std::cout << "totalSize: " << totalSize << "; blockSize: " << blockSize << std::endl;
    throw SGPP::base::operation_exception(
        "totalSize must be divisible by blockSize without remainder, but it is not!");
  }

  // do all further calculations with complete blocks
  size_t blockCount = totalSize / blockSize;

  size_t blockSegmentSize = blockCount / segmentCount;
  size_t remainder = blockCount - blockSegmentSize * segmentCount;
  size_t blockSegmentOffset = 0;

  if (segmentNumber < remainder) {
    blockSegmentSize++;
    blockSegmentOffset = blockSegmentSize * segmentNumber;
  } else {
    blockSegmentOffset =
        remainder * (blockSegmentSize + 1) + (segmentNumber - remainder) * blockSegmentSize;
  }

  *segmentStart = start + blockSegmentOffset* blockSize;
  *segmentEnd = *segmentStart + blockSegmentSize* blockSize;
}

void OperationMultiEvalModMaskStreaming::getOpenMPPartitionSegment(size_t start, size_t end,
                                                                   size_t* segmentStart,
                                                                   size_t* segmentEnd,
                                                                   size_t blocksize) {
  size_t threadCount = omp_get_num_threads();
  size_t myThreadNum = omp_get_thread_num();
  getPartitionSegment(start, end, threadCount, myThreadNum, segmentStart, segmentEnd, blocksize);
}

size_t OperationMultiEvalModMaskStreaming::getChunkGridPoints() {
  // not used by the MIC-implementation
  return 12;
}
size_t OperationMultiEvalModMaskStreaming::getChunkDataPoints() {
#if defined(__MIC__) || defined(__AVX512F__)
  return STREAMING_MODLINEAR_MIC_AVX512_UNROLLING_WIDTH;
#else
  return 24;  // must be divisible by 24
#endif
}

void OperationMultiEvalModMaskStreaming::mult(SGPP::base::DataVector& alpha,
                                              SGPP::base::DataVector& result) {
  this->myTimer_.start();

  size_t originalSize = result.getSize();

  result.resize(this->preparedDataset.getNcols());

  result.setAll(0.0);

#pragma omp parallel
  {
    size_t start;
    size_t end;
    getOpenMPPartitionSegment(0, this->preparedDataset.getNcols(), &start, &end,
                              getChunkDataPoints());

    this->multImpl(this->level, this->index, this->mask, this->offset, &this->preparedDataset,
                   alpha, result, 0, alpha.getSize(), start, end);
  }
  result.resize(originalSize);
  this->duration = this->myTimer_.stop();
}

void OperationMultiEvalModMaskStreaming::multTranspose(SGPP::base::DataVector& source,
                                                       SGPP::base::DataVector& result) {
  this->myTimer_.start();

  size_t originalSize = source.getSize();

  source.resize(this->preparedDataset.getNcols());

  // set padding area to zero
  for (size_t i = originalSize; i < this->preparedDataset.getNcols(); i++) {
    source[i] = 0.0;
  }

  result.setAll(0.0);

#pragma omp parallel
  {
    size_t start;
    size_t end;

    getOpenMPPartitionSegment(0, this->storage->size(), &start, &end, 1);

    this->multTransposeImpl(this->level, this->index, this->mask, this->offset,
                            &this->preparedDataset, source, result, start, end, 0,
                            this->preparedDataset.getNcols());
  }
  source.resize(originalSize);
  this->duration = this->myTimer_.stop();
}

size_t OperationMultiEvalModMaskStreaming::padDataset(SGPP::base::DataMatrix& dataset) {
  size_t vecWidth = this->getChunkDataPoints();

  // Assure that data has a even number of instances -> padding might be needed
  size_t remainder = dataset.getNrows() % vecWidth;
  size_t loopCount = vecWidth - remainder;

  if (loopCount != vecWidth) {
    SGPP::base::DataVector lastRow(dataset.getNcols());
    size_t oldSize = dataset.getNrows();
    dataset.getRow(dataset.getNrows() - 1, lastRow);
    dataset.resize(dataset.getNrows() + loopCount);

    for (size_t i = 0; i < loopCount; i++) {
      dataset.setRow(oldSize + i, lastRow);
    }
  }

  return dataset.getNrows();
}

float_t OperationMultiEvalModMaskStreaming::getDuration() { return this->duration; }

void OperationMultiEvalModMaskStreaming::prepare() { this->recalculateLevelIndexMask(); }

void OperationMultiEvalModMaskStreaming::recalculateLevelIndexMask() {
  // TODO(pfandedd): does the padding work? test
  //    uint32_t localWorkSize = 24;
  size_t localWorkSize = this->getChunkGridPoints();

  size_t remainder = this->storage->size() % localWorkSize;
  size_t padding = 0;

  if (remainder != 0) {
    padding = localWorkSize - remainder;
  }

  size_t gridSize = this->storage->size() + padding;
  size_t dims = this->storage->dim();

  SGPP::base::HashGridIndex::level_type curLevel;
  SGPP::base::HashGridIndex::index_type curIndex;

  // TODO(pfandedd): update the other kernels with this style

  this->level = std::vector<double>(gridSize * dims);
  this->index = std::vector<double>(gridSize * dims);
  this->mask = std::vector<double>(gridSize * dims);
  this->offset = std::vector<double>(gridSize * dims);

  for (size_t i = 0; i < this->storage->size(); i++) {
    for (size_t dim = 0; dim < dims; dim++) {
      storage->get(i)->get(dim, curLevel, curIndex);

      if (curLevel == 1) {
        this->level[i * dims + dim] = 0.0;
        this->index[i * dims + dim] = 0.0;

        uint64_t intmask = 0x0000000000000000;
        this->mask[i * dims + dim] = *reinterpret_cast<double*>(&intmask);

        this->offset[i * dims + dim] = 1.0;
      } else if (curIndex == 1) {
        this->level[i * dims + dim] =
            static_cast<double>(-1.0) * static_cast<double>(1 << curLevel);
        this->index[i * dims + dim] = 0.0;

        uint64_t intmask = 0x0000000000000000;
        this->mask[i * dims + dim] = *reinterpret_cast<double*>(&intmask);

        this->offset[i * dims + dim] = 2.0;
      } else if (curIndex ==
                 static_cast<SGPP::base::HashGridIndex::level_type>(((1 << curLevel) - 1))) {
        this->level[i * dims + dim] = static_cast<double>(1 << curLevel);
        this->index[i * dims + dim] = static_cast<double>(curIndex);

        uint64_t intmask = 0x0000000000000000;
        this->mask[i * dims + dim] = *reinterpret_cast<double*>(&intmask);

        this->offset[i * dims + dim] = 1.0;
      } else {
        this->level[i * dims + dim] = static_cast<double>(1 << curLevel);
        this->index[i * dims + dim] = static_cast<double>(curIndex);

        uint64_t intmask = 0x8000000000000000;
        this->mask[i * dims + dim] = *reinterpret_cast<double*>(&intmask);

        this->offset[i * dims + dim] = 1.0;
      }
    }
  }

  for (size_t i = this->storage->size(); i < gridSize; i++) {
    for (size_t dim = 0; dim < dims; dim++) {
      this->level[i * dims + dim] = 0;
      this->index[i * dims + dim] = 0;

      uint64_t intmask = 0x0000000000000000;
      this->mask[i * dims + dim] = *reinterpret_cast<double*>(&intmask);

      this->offset[i * dims + dim] = 1.0;
    }
  }
}
}  // namespace datadriven
}  // namespace SGPP

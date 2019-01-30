// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultiEvalStreamingBSpline/OperationMultiEvalStreamingBSpline.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

OperationMultiEvalStreamingBSpline::OperationMultiEvalStreamingBSpline(base::Grid& grid,
                                                                       base::DataMatrix& dataset)
    : OperationMultipleEval(grid, dataset),
      preparedDataset(dataset),
      myTimer_(sgpp::base::SGppStopwatch()),
      duration(-1.0) {
  this->storage = &grid.getStorage();

  this->padDataset(this->preparedDataset);
  this->preparedDataset.transpose();

  // create the kernel specific data structures for the current grid
  this->prepare();
}

OperationMultiEvalStreamingBSpline::~OperationMultiEvalStreamingBSpline() {
  if (this->level_ != nullptr) delete this->level_;

  if (this->index_ != nullptr) delete this->index_;
}

void OperationMultiEvalStreamingBSpline::getPartitionSegment(size_t start, size_t end,
                                                             size_t segmentCount,
                                                             size_t segmentNumber,
                                                             size_t* segmentStart,
                                                             size_t* segmentEnd, size_t blockSize) {
  size_t totalSize = end - start;

  // check for valid input
  if (blockSize == 0) {
    throw sgpp::base::operation_exception("blockSize must not be zero!");
  }

  if (totalSize % blockSize != 0) {
    // std::cout << "totalSize: " << totalSize << "; blockSize: " << blockSize << std::endl;
    throw sgpp::base::operation_exception(
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

  *segmentStart = start + blockSegmentOffset * blockSize;
  *segmentEnd = *segmentStart + blockSegmentSize * blockSize;
}

void OperationMultiEvalStreamingBSpline::getOpenMPPartitionSegment(size_t start, size_t end,
                                                                   size_t* segmentStart,
                                                                   size_t* segmentEnd,
                                                                   size_t blocksize) {
  size_t threadCount = omp_get_num_threads();
  size_t myThreadNum = omp_get_thread_num();
  getPartitionSegment(start, end, threadCount, myThreadNum, segmentStart, segmentEnd, blocksize);
}

size_t OperationMultiEvalStreamingBSpline::getChunkGridPoints() {
  // not used by the MIC-implementation
  return 12;
}
size_t OperationMultiEvalStreamingBSpline::getChunkDataPoints() {
  //#else
  return 8;  // must be divisible by 8
  //#endif
}

void OperationMultiEvalStreamingBSpline::mult(sgpp::base::DataVector& alpha,
                                              sgpp::base::DataVector& result) {
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

    this->multImpl(level_, index_, &this->preparedDataset, alpha, result, 0, alpha.getSize(), start,
                   end);
  }
  result.resize(originalSize);
  this->duration = this->myTimer_.stop();
}

void OperationMultiEvalStreamingBSpline::multTranspose(sgpp::base::DataVector& source,
                                                       sgpp::base::DataVector& result) {
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

    getOpenMPPartitionSegment(0, this->storage->getSize(), &start, &end, 1);

    this->multTransposeImpl(this->level_, this->index_, &this->preparedDataset, source, result,
                            start, end, 0, this->preparedDataset.getNcols());
  }
  source.resize(originalSize);
  this->duration = this->myTimer_.stop();
}

void OperationMultiEvalStreamingBSpline::recalculateLevelAndIndex() {
  if (this->level_ != nullptr) delete this->level_;

  if (this->index_ != nullptr) delete this->index_;

  this->level_ =
      new sgpp::base::DataMatrix(this->storage->getSize(), this->storage->getDimension());
  this->index_ =
      new sgpp::base::DataMatrix(this->storage->getSize(), this->storage->getDimension());

  this->storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
}

size_t OperationMultiEvalStreamingBSpline::padDataset(sgpp::base::DataMatrix& dataset) {
  size_t vecWidth = this->getChunkDataPoints();

  // Assure that data has a even number of instances -> padding might be needed
  size_t remainder = dataset.getNrows() % vecWidth;
  size_t loopCount = vecWidth - remainder;

  if (loopCount != vecWidth) {
    sgpp::base::DataVector lastRow(dataset.getNcols());
    size_t oldSize = dataset.getNrows();
    dataset.getRow(dataset.getNrows() - 1, lastRow);
    dataset.resize(dataset.getNrows() + loopCount);

    for (size_t i = 0; i < loopCount; i++) {
      dataset.setRow(oldSize + i, lastRow);
    }
  }

  return dataset.getNrows();
}

double OperationMultiEvalStreamingBSpline::getDuration() { return this->duration; }

void OperationMultiEvalStreamingBSpline::prepare() { this->recalculateLevelAndIndex(); }
}  // namespace datadriven
}  // namespace sgpp

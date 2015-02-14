// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultiEvalStreaming/OperationMultiEvalStreaming.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

#undef __SSE3__
#undef __AVX__

OperationMultiEvalStreaming::OperationMultiEvalStreaming(base::Grid &grid,
        base::DataMatrix &dataset) :
    OperationMultipleEval(grid, dataset), preparedDataset(dataset), myTimer_(SGPP::base::SGppStopwatch()), duration(-1.0) {
    this->storage = grid.getStorage();
    this->padDataset(this->preparedDataset);
    this->preparedDataset.transpose();
}

OperationMultiEvalStreaming::~OperationMultiEvalStreaming() {
    if (this->level_ != nullptr)
        delete this->level_;

    if (this->index_ != nullptr)
        delete this->index_;
}

void OperationMultiEvalStreaming::getPartitionSegment(size_t start, size_t end, size_t segmentCount,
        size_t segmentNumber, size_t* segmentStart, size_t* segmentEnd, size_t blockSize) {
    size_t totalSize = end - start;

    // check for valid input
    if (blockSize == 0) {
        throw SGPP::base::operation_exception("blockSize must not be zero!");
    }

    if (totalSize % blockSize != 0) {
        //std::cout << "totalSize: " << totalSize << "; blockSize: " << blockSize << std::endl;
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
        blockSegmentOffset = remainder * (blockSegmentSize + 1) + (segmentNumber - remainder) * blockSegmentSize;
    }

    *segmentStart = start + blockSegmentOffset * blockSize;
    *segmentEnd = *segmentStart + blockSegmentSize * blockSize;
}

void OperationMultiEvalStreaming::getOpenMPPartitionSegment(size_t start, size_t end,
        size_t* segmentStart, size_t* segmentEnd, size_t blocksize) {
    size_t threadCount = omp_get_num_threads();
    size_t myThreadNum = omp_get_thread_num();
    getPartitionSegment(start, end, threadCount, myThreadNum, segmentStart, segmentEnd, blocksize);
}

size_t OperationMultiEvalStreaming::getChunkGridPoints() {
    return 12;
}
size_t OperationMultiEvalStreaming::getChunkDataPoints() {
    return 24; // must be divisible by 24
}

void OperationMultiEvalStreaming::mult(SGPP::base::DataVector &alpha, SGPP::base::DataVector &result) {
    this->myTimer_.start();

    this->recalculateLevelAndIndex();

    size_t originalSize = result.getSize();

    result.resize(this->preparedDataset.getNcols());

    result.setAll(0.0);

    #pragma omp parallel
    {
        size_t start;
        size_t end;
        getOpenMPPartitionSegment(0, this->preparedDataset.getNcols(), &start, &end, getChunkDataPoints());
        this->multImpl(level_, index_, &this->preparedDataset, alpha, result, 0, alpha.getSize(), start, end);
    }
    result.resize(originalSize);
    this->duration = this->myTimer_.stop();
}

void OperationMultiEvalStreaming::multTranspose(SGPP::base::DataVector &source,
        SGPP::base::DataVector &result) {
    this->myTimer_.start();
    this->recalculateLevelAndIndex();

    size_t originalSize = source.getSize();

    source.resize(this->preparedDataset.getNcols());
    //set padding area to zero
    for (size_t i = originalSize; i < this->preparedDataset.getNcols(); i++) {
        source[i] = 0.0;
    }

    result.setAll(0.0);

    #pragma omp parallel
    {
        size_t start;
        size_t end;

        getOpenMPPartitionSegment(0, this->storage->size(), &start, &end, 1);
        this->multTransposeImpl(this->level_, this->index_, &this->preparedDataset, source, result, start, end, 0,
                                this->preparedDataset.getNcols());

    }
    source.resize(originalSize);
    this->duration = this->myTimer_.stop();
}

void OperationMultiEvalStreaming::recalculateLevelAndIndex() {
    if (this->level_ != nullptr)
        delete this->level_;

    if (this->index_ != nullptr)
        delete this->index_;

    this->level_ = new SGPP::base::DataMatrix(this->storage->size(), this->storage->dim());
    this->index_ = new SGPP::base::DataMatrix(this->storage->size(), this->storage->dim());

    this->storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
}

size_t OperationMultiEvalStreaming::padDataset(SGPP::base::DataMatrix &dataset) {

    size_t vecWidth = 24;

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

float_t OperationMultiEvalStreaming::getDuration() {
    return this->duration;
}

}
}

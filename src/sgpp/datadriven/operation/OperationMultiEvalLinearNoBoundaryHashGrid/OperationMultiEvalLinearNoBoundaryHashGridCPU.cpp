#include "datadriven/operation/OperationMultiEvalLinearNoBoundaryHashGrid/OperationMultiEvalLinearNoBoundaryHashGridCPU.hpp"

namespace sg {
namespace datadriven {

#undef __SSE3__
#undef __AVX__

OperationMultiEvalLinearNoBoundaryHashGridCPU::OperationMultiEvalLinearNoBoundaryHashGridCPU(base::Grid &grid,
		base::DataMatrix &dataset) :
		OperationMultipleEval(grid, dataset), myTimer_(sg::base::SGppStopwatch()), duration(-1.0) {
	this->storage = grid.getStorage();
}

OperationMultiEvalLinearNoBoundaryHashGridCPU::~OperationMultiEvalLinearNoBoundaryHashGridCPU() {
	if (this->level_ != nullptr)
		delete this->level_;

	if (this->index_ != nullptr)
		delete this->index_;
}

void OperationMultiEvalLinearNoBoundaryHashGridCPU::getPartitionSegment(size_t start, size_t end, size_t segmentCount,
		size_t segmentNumber, size_t* segmentStart, size_t* segmentEnd, size_t blockSize) {
	size_t totalSize = end - start;

	// check for valid input
	if (blockSize == 0) {
		throw sg::base::operation_exception("blockSize must not be zero!");
	}

	if (totalSize % blockSize != 0) {
		//std::cout << "totalSize: " << totalSize << "; blockSize: " << blockSize << std::endl;
		throw sg::base::operation_exception(
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

void OperationMultiEvalLinearNoBoundaryHashGridCPU::getOpenMPPartitionSegment(size_t start, size_t end,
		size_t* segmentStart, size_t* segmentEnd, size_t blocksize) {
	size_t threadCount = 1;
	size_t myThreadNum = 0;
#ifdef _OPENMP
	threadCount = omp_get_num_threads();
	myThreadNum = omp_get_thread_num();
#endif
	getPartitionSegment(start, end, threadCount, myThreadNum, segmentStart, segmentEnd, blocksize);
}

size_t OperationMultiEvalLinearNoBoundaryHashGridCPU::getChunkGridPoints() {
	return 12;
}
size_t OperationMultiEvalLinearNoBoundaryHashGridCPU::getChunkDataPoints() {
	return 24; // must be divisible by 24
}

void OperationMultiEvalLinearNoBoundaryHashGridCPU::mult(sg::base::DataVector &alpha, sg::base::DataVector &result) {
	this->myTimer_.start();

	this->recalculateLevelAndIndex();

	size_t originalSize = result.getSize();

	sg::base::DataMatrix tmpDataSet(dataset);
	size_t paddedSize = this->padDataset(tmpDataSet);
	tmpDataSet.transpose();

	result.resize(paddedSize);

	result.setAll(0.0);

#pragma omp parallel
	{
//			m_kernel.mult(level_, index_, mask_, offset_, dataset_, alpha,
//					result, 0, alpha.getSize(), m_datasetFrom, m_datasetTo);
		size_t start;
		size_t end;
		getOpenMPPartitionSegment(0, tmpDataSet.getNcols(), &start, &end, getChunkDataPoints());
		this->multImpl(level_, index_, &tmpDataSet, alpha, result, 0, alpha.getSize(), start, end);
	}
	result.resize(originalSize);
	//return this->myTimer_.stop();
	this->duration = this->myTimer_.stop();
}

void OperationMultiEvalLinearNoBoundaryHashGridCPU::multTranspose(sg::base::DataVector &source,
		sg::base::DataVector &result) {
	this->myTimer_.start();
	this->recalculateLevelAndIndex();

	size_t originalSize = source.getSize();

	sg::base::DataMatrix tmpDataSet(dataset);
	size_t paddedSize = this->padDataset(tmpDataSet);
	tmpDataSet.transpose();

	source.resize(paddedSize);
	//set padding area to zero
	for (size_t i = originalSize; i < paddedSize; i++) {
		source[i] = 0.0;
	}

	result.setAll(0.0);

#pragma omp parallel
	{
		size_t start;
		size_t end;

		getOpenMPPartitionSegment(0, this->storage->size(), &start, &end, 1);
		this->multTransposeImpl(this->level_, this->index_, &tmpDataSet, source, result, start, end, 0,
				tmpDataSet.getNcols());

	}
	source.resize(originalSize);
	this->duration = this->myTimer_.stop();
	//return this->myTimer_.stop();
}

void OperationMultiEvalLinearNoBoundaryHashGridCPU::recalculateLevelAndIndex() {
	if (this->level_ != nullptr)
		delete this->level_;

	if (this->index_ != nullptr)
		delete this->index_;

	this->level_ = new sg::base::DataMatrix(this->storage->size(), this->storage->dim());
	this->index_ = new sg::base::DataMatrix(this->storage->size(), this->storage->dim());

	this->storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
}

size_t OperationMultiEvalLinearNoBoundaryHashGridCPU::padDataset(sg::base::DataMatrix &dataset) {

	size_t vecWidth = 24;

	// Assure that data has a even number of instances -> padding might be needed
	size_t remainder = dataset.getNrows() % vecWidth;
	size_t loopCount = vecWidth - remainder;

	if (loopCount != vecWidth) {
		sg::base::DataVector lastRow(dataset.getNcols());
		size_t oldSize = dataset.getNrows();
		dataset.getRow(dataset.getNrows() - 1, lastRow);
		dataset.resize(dataset.getNrows() + loopCount);

		for (size_t i = 0; i < loopCount; i++) {
			dataset.setRow(oldSize + i, lastRow);
		}
	}

	return dataset.getNrows();
}

double OperationMultiEvalLinearNoBoundaryHashGridCPU::getDuration() {
	return this->duration;
}

}
}

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OperationMultiEvalStreamingOCL.hpp"

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

OperationMultiEvalStreamingOCL::OperationMultiEvalStreamingOCL(base::Grid& grid,
		base::DataMatrix& dataset) :
		OperationMultipleEval(grid, dataset), preparedDataset(dataset), myTimer_(
		SGPP::base::SGppStopwatch()), duration(-1.0) { //, mask_(null_ptr), offset_(null_ptr)
	this->manager.initializePlattform();
	this->kernel = new OCLKernelImpl(this->manager);

	this->storage = grid.getStorage();
	this->padDataset(this->preparedDataset);
	this->preparedDataset.transpose();
}

OperationMultiEvalStreamingOCL::~OperationMultiEvalStreamingOCL() {
	if (this->level_ != nullptr)
		delete this->level_;

	if (this->index_ != nullptr)
		delete this->index_;
}

//size_t OperationMultiEvalStreamingOCL::getChunkGridPoints() {
//	return 12;
//}
//size_t OperationMultiEvalStreamingOCL::getChunkDataPoints() {
//	return 24; // must be divisible by 24
//}

void OperationMultiEvalStreamingOCL::mult(SGPP::base::DataVector& alpha,
SGPP::base::DataVector& result) {
	this->myTimer_.start();

	//TODO: do so only if necessary, also in the other kernels
	this->recalculateLevelAndIndex();

	size_t originalSize = result.getSize();

	result.resize(this->preparedDataset.getNcols());

	result.setAll(0.0);

	size_t gridFrom = 0;
	size_t gridTo = grid.getStorage()->size();
	size_t datasetFrom = 0;
	size_t datasetTo = this->preparedDataset.getNcols();

	this->kernel->mult(level_, index_, &this->preparedDataset, alpha, result,
			gridFrom, gridTo, datasetFrom, datasetTo);

	result.resize(originalSize);
	this->duration = this->myTimer_.stop();
}

void OperationMultiEvalStreamingOCL::multTranspose(
SGPP::base::DataVector& source,
SGPP::base::DataVector& result) {
	this->myTimer_.start();
	this->recalculateLevelAndIndex();

	size_t originalSize = source.getSize();
	size_t gridOriginalSize = result.getSize();
	result.resize(this->level_->getNrows());

	source.resize(this->preparedDataset.getNcols());

	//set padding area to zero
	for (size_t i = originalSize; i < this->preparedDataset.getNcols(); i++) {
		source[i] = 0.0;
	}

	result.setAll(0.0);

	size_t gridFrom = 0;
	//size_t gridTo = grid.getStorage()->size();
	size_t gridTo = this->level_->getNrows();
	size_t datasetFrom = 0;
	size_t datasetTo = this->preparedDataset.getNcols();

	this->kernel->multTranspose(this->level_, this->index_,
			&this->preparedDataset, source, result, gridFrom, gridTo,
			datasetFrom, datasetTo);

	source.resize(originalSize);
	result.resize(gridOriginalSize);
	this->duration = this->myTimer_.stop();
}

void OperationMultiEvalStreamingOCL::recalculateLevelAndIndex() {
	if (this->level_ != nullptr)
		delete this->level_;

	if (this->index_ != nullptr)
		delete this->index_;

	uint32_t localWorkSize = this->manager.getOCLLocalSize();

	size_t remainder = this->storage->size() % localWorkSize;
	size_t padding = localWorkSize - remainder;
	size_t paddedSize = this->storage->size() + padding;

	this->level_ = new SGPP::base::DataMatrix(paddedSize,
			this->storage->dim());
	this->index_ = new SGPP::base::DataMatrix(paddedSize,
			this->storage->dim());

	this->storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));


	for (size_t i = this->storage->size(); i < paddedSize; i++) {
		for (size_t j = 0; j < this->storage->dim(); j++) {
			this->level_->set(i, j, 1.0);
			this->index_->set(i, j, 1.0);
		}
	}
	//TODO: one call per CG iteration -> very expensive!
	this->kernel->resetKernel();
}

//TODO: fix for OpenCL
size_t OperationMultiEvalStreamingOCL::padDataset(
SGPP::base::DataMatrix& dataset) {

	size_t vecWidth = this->manager.getOCLLocalSize();

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

float_t OperationMultiEvalStreamingOCL::getDuration() {
	return this->duration;
}

}
}

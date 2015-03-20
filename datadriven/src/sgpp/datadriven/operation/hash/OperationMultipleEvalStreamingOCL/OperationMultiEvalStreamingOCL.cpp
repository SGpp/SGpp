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
		OperationMultipleEval(grid, dataset), preparedDataset(dataset), myTimer(
		SGPP::base::SGppStopwatch()), duration(-1.0) {
	this->manager.initializePlattform();
	this->dims = dataset.getNcols(); //be aware of transpose!
	this->kernel = new OCLKernelImpl<double>(dims, this->manager);

	this->storage = grid.getStorage();
	this->padDataset(this->preparedDataset);
	this->preparedDataset.transpose();
	this->datasetSize = this->preparedDataset.getNcols();

	std::cout << "dims: " << this->dims << std::endl;
	std::cout << "padded instances: " << this->datasetSize << std::endl;

	this->kernelDataset = new double[this->preparedDataset.getNrows()
			* this->preparedDataset.getNcols()];
	for (size_t i = 0; i < this->preparedDataset.getSize(); i++) {
		this->kernelDataset[i] = (double) this->preparedDataset[i];
	}

	//create the kernel specific data structures
	this->prepare();
}

OperationMultiEvalStreamingOCL::~OperationMultiEvalStreamingOCL() {
	if (this->level != nullptr) {
		delete this->level;
	}

	if (this->index != nullptr) {
		delete this->index;
	}

	if (this->kernelDataset != nullptr) {
		delete this->kernelDataset;
	}
}

void OperationMultiEvalStreamingOCL::mult(SGPP::base::DataVector& alpha,
SGPP::base::DataVector& result) {
	this->myTimer.start();

	size_t originalSize = result.getSize();

	result.resize(this->preparedDataset.getNcols());
	result.setAll(0.0);

	size_t originalAlphaSize = alpha.getSize();
	alpha.resize(this->gridSize);
	for (size_t i = originalAlphaSize; i < alpha.getSize(); i++) {
		alpha[i] = 0.0;
	}

	size_t gridFrom = 0;
	size_t gridTo = this->gridSize;
	size_t datasetFrom = 0;
	size_t datasetTo = this->datasetSize;

	this->kernel->mult(this->level, this->index, this->gridSize,
			this->kernelDataset, this->datasetSize, alpha, result, gridFrom,
			gridTo, datasetFrom, datasetTo);

	result.resize(originalSize);
	alpha.resize(originalAlphaSize);
	this->duration = this->myTimer.stop();
}

void OperationMultiEvalStreamingOCL::multTranspose(
SGPP::base::DataVector& source,
SGPP::base::DataVector& result) {
	this->myTimer.start();

	size_t originalSize = source.getSize();
	size_t gridOriginalSize = result.getSize();
	result.resize(this->gridSize);

	source.resize(this->preparedDataset.getNcols());

	//set padding area to zero
	for (size_t i = originalSize; i < this->preparedDataset.getNcols(); i++) {
		source[i] = 0.0;
	}

	result.setAll(0.0);

	size_t gridFrom = 0;
	//size_t gridTo = grid.getStorage()->size();
	size_t gridTo = this->gridSize;
	size_t datasetFrom = 0;
	size_t datasetTo = this->preparedDataset.getNcols();

	this->kernel->multTranspose(this->level, this->index, this->gridSize,
			this->kernelDataset, this->preparedDataset.getNcols(), source,
			result, gridFrom, gridTo, datasetFrom, datasetTo);

	source.resize(originalSize);
	result.resize(gridOriginalSize);
	this->duration = this->myTimer.stop();
}

void OperationMultiEvalStreamingOCL::recalculateLevelAndIndex() {
	if (this->level != nullptr)
		delete this->level;

	if (this->index != nullptr)
		delete this->index;

	uint32_t localWorkSize = this->manager.getOCLLocalSize();

	size_t remainder = this->storage->size() % localWorkSize;
	size_t padding = 0;
	if (remainder != 0) {
		padding = localWorkSize - remainder;
	}
	this->gridSize = this->storage->size() + padding;

	SGPP::base::DataMatrix *levelMatrix = new SGPP::base::DataMatrix(
			this->gridSize, this->dims);
	SGPP::base::DataMatrix *indexMatrix = new SGPP::base::DataMatrix(
			this->gridSize, this->dims);

	this->storage->getLevelIndexArraysForEval(*levelMatrix, *indexMatrix);

	for (size_t i = this->storage->size(); i < this->gridSize; i++) {
		for (size_t j = 0; j < this->storage->dim(); j++) {
			levelMatrix->set(i, j, 1.0);
			indexMatrix->set(i, j, 1.0);
		}
	}

	this->level = new double[this->gridSize * this->dims];
	this->index = new double[this->gridSize * this->dims];

	for (size_t i = 0; i < this->gridSize * this->dims; i++) {
		this->level[i] = (*levelMatrix)[i];
		this->index[i] = (*indexMatrix)[i];
	}

	delete levelMatrix;
	delete indexMatrix;

}

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

void OperationMultiEvalStreamingOCL::prepare() {
	this->recalculateLevelAndIndex();

	this->kernel->resetKernel();

	std::cout << "gridSize: " << this->gridSize << std::endl;
}

}
}

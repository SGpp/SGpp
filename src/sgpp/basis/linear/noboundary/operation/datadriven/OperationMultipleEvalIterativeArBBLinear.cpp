/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeArBBLinear.hpp"
#include "exception/operation_exception.hpp"
using namespace sg::base;

namespace sg
{
namespace parallel
{

OperationMultipleEvalIterativeArBBLinear::OperationMultipleEvalIterativeArBBLinear(GridStorage* storage, DataMatrix* dataset) : OperationMultipleEvalVectorized(dataset)
{
	this->storage = storage;

	this->level_ = new DataMatrix(storage->size(), storage->dim());
	this->index_ = new DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myTimer = new SGppStopwatch();
	myArBBKernels = new ArBBKernels();
}

OperationMultipleEvalIterativeArBBLinear::~OperationMultipleEvalIterativeArBBLinear()
{
	delete myTimer;
	delete myArBBKernels;
}

void OperationMultipleEvalIterativeArBBLinear::rebuildLevelAndIndex()
{
	delete this->level_;
	delete this->index_;

	this->level_ = new DataMatrix(storage->size(), storage->dim());
	this->index_ = new DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myArBBKernels->resetKernels();
}

double OperationMultipleEvalIterativeArBBLinear::multVectorized(DataVector& alpha, DataVector& result)
{
	size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    result.setAll(0.0);

    double* ptrAlpha = alpha.getPointer();
    double* ptrData =  this->dataset_->getPointer();
    double* ptrLevel = this->level_->getPointer();
    double* ptrIndex = this->index_->getPointer();
    double* ptrResult = result.getPointer();

    if (this->dataset_->getNrows() % 16 != 0 || result_size !=  this->dataset_->getNrows())
    {
    	throw operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
    }

    double time = myArBBKernels->multArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);

	return time;
}

double OperationMultipleEvalIterativeArBBLinear::multTransposeVectorized(DataVector& source, DataVector& result)
{
	size_t soruceSize = source.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    result.setAll(0.0);

    double* ptrSource = source.getPointer();
    double* ptrData = this->dataset_->getPointer();
    double* ptrLevel = this->level_->getPointer();
    double* ptrIndex = this->index_->getPointer();
    double* ptrGlobalResult = result.getPointer();

    if (this->dataset_->getNrows() % 16 != 0 || soruceSize != this->dataset_->getNrows())
    {
    	throw operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    double time = myArBBKernels->multTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, soruceSize, storageSize, dims);

   	return time;
}

}
}

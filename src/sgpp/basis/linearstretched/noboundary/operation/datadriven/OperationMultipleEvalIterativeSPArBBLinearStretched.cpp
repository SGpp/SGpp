/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPArBBLinearStretched.hpp"
#include "exception/operation_exception.hpp"

namespace sg
{

OperationMultipleEvalIterativeSPArBBLinearStretched::OperationMultipleEvalIterativeSPArBBLinearStretched(GridStorage* storage, DataMatrixSP* dataset) : OperationMultipleEvalVectorizedSP(dataset)
{
	this->storage = storage;

	this->level_ = new DataMatrixSP(storage->size(), storage->dim());
	this->index_ = new DataMatrixSP(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myTimer = new SGppStopwatch();
	myArBBKernels = new ArBBKernels();
}

OperationMultipleEvalIterativeSPArBBLinearStretched::~OperationMultipleEvalIterativeSPArBBLinearStretched()
{
	delete myTimer;
	delete myArBBKernels;
}

void OperationMultipleEvalIterativeSPArBBLinearStretched::rebuildLevelAndIndex()
{
	delete this->level_;
	delete this->index_;

	this->level_ = new DataMatrixSP(storage->size(), storage->dim());
	this->index_ = new DataMatrixSP(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myArBBKernels->resetKernels();
}

double OperationMultipleEvalIterativeSPArBBLinearStretched::multTransposeVectorized(DataVectorSP& source, DataVectorSP& result)
{
	size_t source_size = source.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    result.setAll(0.0f);

    float* ptrSource = source.getPointer();
    float* ptrData = this->dataset_->getPointer();
    float* ptrLevel = this->level_->getPointer();
    float* ptrIndex = this->index_->getPointer();
    float* ptrGlobalResult = result.getPointer();

    if (this->dataset_->getNrows() % 16 != 0 || source_size != this->dataset_->getNrows())
    {
    	throw operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
    }

    double time = myArBBKernels->multTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);

	return time;
}

double OperationMultipleEvalIterativeSPArBBLinearStretched::multVectorized(DataVectorSP& alpha, DataVectorSP& result)
{
	size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    result.setAll(0.0f);

    float* ptrAlpha = alpha.getPointer();
    float* ptrData = this->dataset_->getPointer();
    float* ptrResult = result.getPointer();
    float* ptrLevel = this->level_->getPointer();
    float* ptrIndex = this->index_->getPointer();

    if (this->dataset_->getNrows() % 16 != 0 || result_size != this->dataset_->getNrows())
    {
    	throw operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    double time = myArBBKernels->multSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);

   	return time;
}

}

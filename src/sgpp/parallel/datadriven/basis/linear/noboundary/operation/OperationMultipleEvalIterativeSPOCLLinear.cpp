/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeSPOCLLinear.hpp"
#include "base/exception/operation_exception.hpp"

namespace sg
{
namespace parallel
{

OperationMultipleEvalIterativeSPOCLLinear::OperationMultipleEvalIterativeSPOCLLinear(
		sg::base::GridStorage* storage, sg::base::DataMatrixSP* dataset) :
	sg::parallel::OperationMultipleEvalVectorizedSP(storage, dataset)
{
	this->level_ = new sg::base::DataMatrixSP(this->storage_->size(), this->storage_->dim());
	this->index_ = new sg::base::DataMatrixSP(this->storage_->size(), this->storage_->dim());

	this->storage_->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myOCLKernels = new OCLKernels();
}

OperationMultipleEvalIterativeSPOCLLinear::~OperationMultipleEvalIterativeSPOCLLinear()
{
	delete myOCLKernels;
}

void OperationMultipleEvalIterativeSPOCLLinear::rebuildLevelAndIndex()
{
	delete this->level_;
	delete this->index_;

	this->level_ = new sg::base::DataMatrixSP(storage_->size(), storage_->dim());
	this->index_ = new sg::base::DataMatrixSP(storage_->size(), storage_->dim());

	storage_->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myOCLKernels->resetKernels();
}

double OperationMultipleEvalIterativeSPOCLLinear::multTransposeVectorized(sg::base::DataVectorSP& source, sg::base::DataVectorSP& result)
{
	size_t source_size = source.getSize();
	size_t dims = storage_->dim();
	size_t storageSize = storage_->size();

    result.setAll(0.0f);

    float* ptrSource = source.getPointer();
    float* ptrData = this->dataset_->getPointer();
    float* ptrLevel = this->level_->getPointer();
    float* ptrIndex = this->index_->getPointer();
    float* ptrGlobalResult = result.getPointer();

    if (this->dataset_->getNrows() % OCL_SGPP_LOCAL_WORKGROUP_SIZE != 0 || source_size != this->dataset_->getNrows())
    {
    	throw sg::base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
    }

	size_t numWGs = storageSize/OCL_SGPP_LOCAL_WORKGROUP_SIZE;
    size_t global = numWGs*OCL_SGPP_LOCAL_WORKGROUP_SIZE;

    double time = 0.0;
    if (global > 0)
    	time = myOCLKernels->multTransSPOCL(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims, global);

	#pragma omp parallel for
	for (size_t j = global; j < storageSize; j++)
	{
		ptrGlobalResult[j] = 0.0f;

		for (size_t i = 0; i < source_size; i++)
		{
			float curSupport = ptrSource[i];

			for (size_t d = 0; d < dims; d++)
			{
				float eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
				float index_calc = eval - (ptrIndex[(j*dims)+d]);
				float abs = (float)fabs(index_calc);
				float last = 1.0f - abs;
				float localSupport = std::max<float>(last, 0.0f);
				curSupport *= localSupport;
			}

			ptrGlobalResult[j] += curSupport;
		}
	}

	return time;
}

double OperationMultipleEvalIterativeSPOCLLinear::multVectorized(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result)
{
	size_t result_size = result.getSize();
	size_t dims = storage_->dim();
	size_t storageSize = storage_->size();

    result.setAll(0.0f);

    float* ptrAlpha = alpha.getPointer();
    float* ptrData = this->dataset_->getPointer();
    float* ptrResult = result.getPointer();
    float* ptrLevel = this->level_->getPointer();
    float* ptrIndex = this->index_->getPointer();

    if (this->dataset_->getNrows() % OCL_SGPP_LOCAL_WORKGROUP_SIZE != 0 || result_size != this->dataset_->getNrows())
    {
    	throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    double time = myOCLKernels->multSPOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims, result_size);

   	return time;
}

}
}

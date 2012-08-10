/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeOCLLinear.hpp"
#include "base/exception/operation_exception.hpp"

namespace sg
{
namespace parallel
{

OperationMultipleEvalIterativeOCLLinear::OperationMultipleEvalIterativeOCLLinear(
		sg::base::GridStorage* storage, sg::base::DataMatrix* dataset,
		int gridFrom, int gridTo, int datasetFrom, int datasetTo) : sg::parallel::OperationMultipleEvalVectorized(dataset)
{
	m_gridFrom = gridFrom;
	m_gridTo = gridTo;
	m_datasetFrom = datasetFrom;
	m_datasetTo = datasetTo;
	adaptDatasetBoundaries();

	this->storage = storage;

	this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
	this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myTimer = new sg::base::SGppStopwatch();
	myOCLKernels = new OCLKernels();
}

OperationMultipleEvalIterativeOCLLinear::~OperationMultipleEvalIterativeOCLLinear()
{
	delete myTimer;
	delete myOCLKernels;
}

void OperationMultipleEvalIterativeOCLLinear::rebuildLevelAndIndex()
{
	delete this->level_;
	delete this->index_;

	this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
	this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myOCLKernels->resetKernels();
}

void OperationMultipleEvalIterativeOCLLinear::updateGridComputeBoundaries(int gridFrom, int gridTo)
{
	m_gridFrom = gridFrom;
	m_gridTo = gridTo;
}

double OperationMultipleEvalIterativeOCLLinear::multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result)
{
	size_t source_size = source.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    result.setAll(0.0);

    double* ptrSource = source.getPointer();
    double* ptrData = this->dataset_->getPointer();
    double* ptrLevel = this->level_->getPointer();
    double* ptrIndex = this->index_->getPointer();
    double* ptrGlobalResult = result.getPointer();

    if (this->dataset_->getNrows() % 128 != 0 || source_size != this->dataset_->getNrows())
    {
    	throw sg::base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
    }

	double time = myOCLKernels->multTransOCL(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims, storageSize);

    // do the rest...
	size_t numWGs = storageSize/OCL_MULT_N_DATAPREFETCH_BLOCKSIZE_DP;
    size_t global = numWGs*OCL_MULT_N_DATAPREFETCH_BLOCKSIZE_DP;

    if (global == 0)
    {
    	global = storageSize;
    }

	#pragma omp parallel for
	for (size_t j = global; j < storageSize; j++)
	{
		ptrGlobalResult[j] = 0.0f;

		for (size_t i = 0; i < source_size; i++)
		{
			double curSupport = ptrSource[i];

			for (size_t d = 0; d < dims; d++)
			{
				double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
				double index_calc = eval - (ptrIndex[(j*dims)+d]);
				double abs = fabs(index_calc);
				double last = 1.0 - abs;
				double localSupport = std::max<double>(last, 0.0);
				curSupport *= localSupport;
			}

			ptrGlobalResult[j] += curSupport;
		}
	}

	return time;
}

double OperationMultipleEvalIterativeOCLLinear::multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    result.setAll(0.0);

    double* ptrAlpha = alpha.getPointer();
    double* ptrData = this->dataset_->getPointer();
    double* ptrResult = result.getPointer();
    double* ptrLevel = this->level_->getPointer();
    double* ptrIndex = this->index_->getPointer();

    if (this->dataset_->getNrows() % 128 != 0 || result_size != this->dataset_->getNrows())
    {
    	throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

	if (m_datasetFrom % 128 != 0 || m_datasetTo % 128 != 0)
	{
		throw sg::base::operation_exception("The borders of the domain segment have to be aligned to multples of 128.");
	}


	double* ptrDataThisProcess = &(ptrData[m_datasetFrom*dims]);
	size_t gpu_partition = m_datasetTo - m_datasetFrom;

	sg::base::DataVector altResult(result.getSize());
	double* ptrAltResult = altResult.getPointer();
	double* ptrResultThisProcess = &(ptrAltResult[m_datasetFrom]);

	std::cout << "from " << m_datasetFrom << " to " << m_datasetTo << "; gpu_partition: "<< gpu_partition <<"; dims: " << dims <<  std::endl;

	for(int i = 0; i< gpu_partition*dims; i++){
		if(ptrData[m_datasetFrom*dims + i] != ptrDataThisProcess[i]){
			std::cout << "data differs[" << i << "]: (orig) " << ptrData[m_datasetFrom*dims + i] << " != " << ptrDataThisProcess[i] << " (partial)" << std::endl;
			throw sg::base::operation_exception("data fail");
		}
	}

	double time = myOCLKernels->multOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims, result_size);

	myOCLKernels->multOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrAltResult, gpu_partition, storageSize, dims, result_size);

	for(int i = 0; i< gpu_partition; i++){
		if(ptrResult[m_datasetFrom + i] != ptrAltResult[m_datasetFrom + i]){
			std::cout << "results differ[" << i << "]: (orig) " << ptrResult[m_datasetFrom + i] << " != " << ptrAltResult[m_datasetFrom + i] << " (partial)" << std::endl;
			throw sg::base::operation_exception("fail");
		}
	}

   	return time;
}

}
}

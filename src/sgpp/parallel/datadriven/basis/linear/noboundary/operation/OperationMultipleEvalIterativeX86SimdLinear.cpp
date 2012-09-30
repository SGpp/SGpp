/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeX86SimdLinear.hpp"
#include "base/exception/operation_exception.hpp"
#include "parallel/tools/PartitioningTool.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinearMult.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinearMultTranspose.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif



namespace sg
{
namespace parallel
{

OperationMultipleEvalIterativeX86SimdLinear::OperationMultipleEvalIterativeX86SimdLinear(
        sg::base::GridStorage* storage, sg::base::DataMatrix* dataset,
		int gridFrom, int gridTo, int datasetFrom, int datasetTo) : sg::parallel::OperationMultipleEvalVectorized(dataset)
{
	m_gridFrom = gridFrom;
	m_gridTo = gridTo;
	m_datasetFrom = datasetFrom;
	m_datasetTo = datasetTo;

#ifdef _OPENMP
#pragma omp parallel
{
    if (omp_get_thread_num() == 0) {
        std::cout << "using " << omp_get_num_threads() << " OpenMP Threads" << std::endl;
    }
}
#endif

	this->storage = storage;

	this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
	this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	myTimer = new sg::base::SGppStopwatch();
}

OperationMultipleEvalIterativeX86SimdLinear::~OperationMultipleEvalIterativeX86SimdLinear()
{
	delete myTimer;
}

void OperationMultipleEvalIterativeX86SimdLinear::rebuildLevelAndIndex()
{
	delete this->level_;
	delete this->index_;

	this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
	this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
}

void OperationMultipleEvalIterativeX86SimdLinear::updateGridComputeBoundaries(int gridFrom, int gridTo)
{
	m_gridFrom = gridFrom;
	m_gridTo = gridTo;
}


double OperationMultipleEvalIterativeX86SimdLinear::multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result)
{
	if (this->dataset_->getNcols() % sg::parallel::X86SimdLinearMultTranspose::getChunkDataPoints() != 0 || source.getSize() != this->dataset_->getNcols())
    {
    	throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

    myTimer->start();
    result.setAll(0.0);

#ifdef _OPENMP
    #pragma omp parallel
	{
#endif
		size_t start;
		size_t end;
		sg::parallel::PartitioningTool::getOpenMPLoopPartitionSegment(m_gridFrom, m_gridTo, &start, &end, 1);

		sg::parallel::X86SimdLinearMultTranspose::multTranspose(level_, index_, dataset_, source, result, start, end, 0, this->dataset_->getNcols());
#ifdef _OPENMP
	}
#endif

	return myTimer->stop();
}

double OperationMultipleEvalIterativeX86SimdLinear::multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	if (this->dataset_->getNcols() % sg::parallel::X86SimdLinearMult::getChunkDataPoints() != 0 || result.getSize() != this->dataset_->getNcols())
    {
    	throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

	myTimer->start();
	result.setAll(0.0);

#ifdef _OPENMP
	#pragma omp parallel
	{
#endif
		size_t start;
		size_t end;
		sg::parallel::PartitioningTool::getOpenMPLoopPartitionSegment(m_datasetFrom, m_datasetTo, &start, &end, sg::parallel::X86SimdLinearMult::getChunkDataPoints());

		if (start % sg::parallel::X86SimdLinearMult::getChunkDataPoints() != 0 || end % sg::parallel::X86SimdLinearMult::getChunkDataPoints() != 0)
        {
			std::cout << "start%CHUNKDATAPOINTS_X86: " << start%sg::parallel::X86SimdLinearMult::getChunkDataPoints() << "; end%CHUNKDATAPOINTS_X86: " << end%sg::parallel::X86SimdLinearMult::getChunkDataPoints() << std::endl;
            throw sg::base::operation_exception("processed vector segment must fit to CHUNKDATAPOINTS_X86!");
        }
		sg::parallel::X86SimdLinearMult::mult(level_, index_, dataset_, alpha, result, start, end);

#ifdef _OPENMP
	}
#endif

	return myTimer->stop();
}

}

}

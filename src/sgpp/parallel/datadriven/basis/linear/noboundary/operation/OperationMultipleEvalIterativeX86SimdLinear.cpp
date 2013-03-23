/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeX86SimdLinear.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinear.hpp"
#include "parallel/tools/PartitioningTool.hpp"

namespace sg
{
namespace parallel
{

OperationMultipleEvalIterativeX86SimdLinear::OperationMultipleEvalIterativeX86SimdLinear(
        sg::base::GridStorage* storage, sg::base::DataMatrix* dataset,
		int gridFrom, int gridTo, int datasetFrom, int datasetTo) :
	sg::parallel::OperationMultipleEvalVectorized(storage, dataset)
{
	m_gridFrom = gridFrom;
	m_gridTo = gridTo;
	m_datasetFrom = datasetFrom;
	m_datasetTo = datasetTo;

	rebuildLevelAndIndex();
}

void OperationMultipleEvalIterativeX86SimdLinear::rebuildLevelAndIndex()
{
	LevelIndexMaskOffsetHelper::rebuild<X86SimdLinear::kernelType>(this);
}

void OperationMultipleEvalIterativeX86SimdLinear::updateGridComputeBoundaries(int gridFrom, int gridTo)
{
	m_gridFrom = gridFrom;
	m_gridTo = gridTo;
}

double OperationMultipleEvalIterativeX86SimdLinear::multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result)
{
	myTimer_->start();
    result.setAll(0.0);

    #pragma omp parallel
	{
		size_t start;
		size_t end;
		PartitioningTool::getOpenMPPartitionSegment(m_gridFrom, m_gridTo, &start, &end, 1);

		X86SimdLinear::multTranspose(level_, index_, NULL, NULL, dataset_, source, result, start, end, 0, this->dataset_->getNcols());
	}

	return myTimer_->stop();
}

double OperationMultipleEvalIterativeX86SimdLinear::multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	myTimer_->start();
	result.setAll(0.0);

	#pragma omp parallel
	{
		size_t start;
		size_t end;
		PartitioningTool::getOpenMPPartitionSegment(m_datasetFrom, m_datasetTo, &start, &end, X86SimdLinear::getChunkDataPoints());

		X86SimdLinear::mult(level_, index_, NULL, NULL, dataset_, alpha, result, 0, alpha.getSize(), start, end);
	}

	return myTimer_->stop();
}

}
}

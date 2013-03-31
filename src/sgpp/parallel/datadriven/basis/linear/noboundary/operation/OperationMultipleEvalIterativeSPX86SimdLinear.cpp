/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeSPX86SimdLinear.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/SPX86SimdLinear.hpp"
#include "parallel/tools/PartitioningTool.hpp"

namespace sg
{
namespace parallel
{

OperationMultipleEvalIterativeSPX86SimdLinear::OperationMultipleEvalIterativeSPX86SimdLinear(
		sg::base::GridStorage* storage, sg::base::DataMatrixSP* dataset,
		int gridFrom, int gridTo, int datasetFrom, int datasetTo) :
	OperationMultipleEvalVectorizedSP(storage, dataset)
{
	m_gridFrom = gridFrom;
	m_gridTo = gridTo;
	m_datasetFrom = datasetFrom;
	m_datasetTo = datasetTo;

	rebuildLevelAndIndex();
}

void OperationMultipleEvalIterativeSPX86SimdLinear::rebuildLevelAndIndex()
{
	LevelIndexMaskOffsetHelperSP::rebuild<SPX86SimdLinear::kernelType, OperationMultipleEvalVectorizedSP>(this);
}

void OperationMultipleEvalIterativeSPX86SimdLinear::updateGridComputeBoundaries(int gridFrom, int gridTo)
{
	m_gridFrom = gridFrom;
	m_gridTo = gridTo;
}

double OperationMultipleEvalIterativeSPX86SimdLinear::multTransposeVectorized(sg::base::DataVectorSP& source, sg::base::DataVectorSP& result)
{
	myTimer_->start();
    result.setAll(0.0f);

	#pragma omp parallel
	{
		size_t start;
		size_t end;
		PartitioningTool::getOpenMPPartitionSegment(m_gridFrom, m_gridTo, &start, &end, 1);

		SPX86SimdLinear::multTranspose(level_, index_, NULL, NULL, dataset_, source, result, start, end, 0, dataset_->getNcols());
	}

	return myTimer_->stop();
}

double OperationMultipleEvalIterativeSPX86SimdLinear::multVectorized(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result)
{
	myTimer_->start();
	result.setAll(0.0f);

	#pragma omp parallel
	{
		size_t start;
		size_t end;
		PartitioningTool::getOpenMPPartitionSegment(m_datasetFrom, m_datasetTo, &start, &end, SPX86SimdLinear::getChunkDataPoints());

		SPX86SimdLinear::mult(level_, index_, NULL, NULL, dataset_, alpha, result, 0, alpha.getSize(), start, end);
	}

	return myTimer_->stop();
}

}
}

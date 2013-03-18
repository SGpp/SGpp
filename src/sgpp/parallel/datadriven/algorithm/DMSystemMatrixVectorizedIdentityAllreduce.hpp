/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef DMSYSTEMMATRIXVECTORIZEDIDENTITYALLREDUCE_H
#define DMSYSTEMMATRIXVECTORIZEDIDENTITYALLREDUCE_H

#include "datadriven/algorithm/DMSystemMatrixBase.hpp"
#include "parallel/tools/MPI/SGppMPITools.hpp"
#include "parallel/tools/PartitioningTool.hpp"
#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/exception/operation_exception.hpp"
#include "base/grid/Grid.hpp"
#include "parallel/tools/TypesParallel.hpp"

namespace sg
{
namespace parallel
{
template<typename MultType, typename MultTransType>
class DMSystemMatrixVectorizedIdentityAllreduce : public sg::datadriven::DMSystemMatrixBase
{
private:
	VectorizationType vecMode_;
	size_t numTrainingInstances_;
	size_t numPatchedTrainingInstances_;

	// which part of the dataset this process handles
	size_t data_size;
	size_t data_offset;
	sg::base::Grid& m_grid;
	sg::base::DataMatrix* level_;
	/// Member to store the sparse grid's indices for better vectorization
	sg::base::DataMatrix* index_;

	// only allocate temporary arrays once
	sg::base::DataVector *tempData;
	sg::base::DataVector *result_tmp;

public:
	DMSystemMatrixVectorizedIdentityAllreduce(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
		: DMSystemMatrixBase(trainData, lambda), vecMode_(vecMode), numTrainingInstances_(0), numPatchedTrainingInstances_(0), m_grid(SparseGrid)
	{
		// handle unsupported vector extensions
		if (this->vecMode_ != X86SIMD && this->vecMode_ != MIC && this->vecMode_ != Hybrid_X86SIMD_MIC && this->vecMode_ != OpenCL && this->vecMode_ != ArBB && this->vecMode_ != Hybrid_X86SIMD_OpenCL)
		{
			throw new sg::base::operation_exception("DMSystemMatrixVectorizedIdentityAllreduce : un-supported vector extension!");
		}

		this->dataset_ = new sg::base::DataMatrix(trainData);
		this->numTrainingInstances_ = this->dataset_->getNrows();
		this->numPatchedTrainingInstances_ = sg::parallel::DMVectorizationPaddingAssistant::padDataset(*(this->dataset_), vecMode_);
		this->tempData = new sg::base::DataVector(this->numPatchedTrainingInstances_);

		if (this->vecMode_ != OpenCL && this->vecMode_ != ArBB && this->vecMode_ != Hybrid_X86SIMD_OpenCL)
		{
			this->dataset_->transpose();
		}

		sg::parallel::PartitioningTool::getMPIPartitionSegment(numPatchedTrainingInstances_, &data_size, &data_offset, sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));

		this->level_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());
		this->index_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());

		m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
		this->result_tmp = new sg::base::DataVector(m_grid.getSize());
	}

	virtual ~DMSystemMatrixVectorizedIdentityAllreduce(){
		delete this->dataset_;

		delete this->level_;
		delete this->index_;

		delete this->result_tmp;
		delete this->tempData;
	}

	virtual void mult(base::DataVector &alpha, base::DataVector &result){
		tempData->setAll(0.0);
		result_tmp->setAll(0.0);
		result.setAll(0.0);

		this->myTimer_->start();
		#pragma omp parallel
		{
			size_t threadChunkStart, threadChunkEnd;
			sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
					data_offset, data_offset + data_size,
					&threadChunkStart, &threadChunkEnd,
					sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));

			MultType::mult(
					level_,
					index_,
					dataset_,
					alpha,
					*tempData,
					0,
					m_grid.getSize(),
					threadChunkStart,
					threadChunkEnd);

			// patch result -> set additional entries zero
			// only done for processes that need this part of the temp data for multTrans
			for (size_t i = std::max<size_t>(this->numTrainingInstances_, threadChunkStart); i < threadChunkEnd; i++)
			{
				tempData->set(i, 0.0f);
			}

#pragma omp single
			{
				// the time measured here does not represent the complete
				// time spent computing mult, it's probably the first threads
				// that finished mult() that enters here while the other
				// threads might still be busy with mult
				double timeMult = this->myTimer_->stop();
				this->computeTimeMult_ += timeMult;
				this->completeTimeMult_ += timeMult;
				this->myTimer_->start();
			}
			// #pragma omp barrier
			// implicit openMP barrier here (after omp single), which is
			// needed, as multTranspose works on the full data part of this
			// process, so threads might work on unfinished results of mult

			sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
					0, m_grid.getSize(),
					&threadChunkStart, &threadChunkEnd, 1);
			MultTransType::multTranspose(
					level_,
					index_,
					dataset_,
					*tempData,
					*result_tmp,
					threadChunkStart,
					threadChunkEnd,
					data_offset,
					data_offset + data_size);
		}
		this->computeTimeMultTrans_ += this->myTimer_->stop();
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(result_tmp->getPointer(), result.getPointer(), result.getSize(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		this->completeTimeMultTrans_ += this->myTimer_->stop();

//		if(sg::parallel::myGlobalMPIComm->getMyRank() == 0){
//			std::cout << "*";
//		}
	}

	virtual void generateb(base::DataVector &classes, base::DataVector &b){
		tempData->setAll(0.0);
		tempData->copyFrom(classes);
		result_tmp->setAll(0.0);

		#pragma omp parallel
		{
			size_t threadChunkStart, threadChunkEnd;
			sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
					0, m_grid.getSize(),
					&threadChunkStart, &threadChunkEnd, 1);
			MultTransType::multTranspose(
						level_,
						index_,
						dataset_,
						*tempData,
						*result_tmp,
						threadChunkStart,
						threadChunkEnd,
						data_offset,
						data_offset + data_size);
		}
		MPI_Allreduce(result_tmp->getPointer(), b.getPointer(), b.getSize(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}

	virtual void rebuildLevelAndIndex(){
		delete this->level_;
		delete this->index_;
		delete this->result_tmp;

		this->level_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());
		this->index_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());

		m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
		result_tmp = new sg::base::DataVector(m_grid.getSize());
	}

};

}
}

#endif // DMSYSTEMMATRIXVECTORIZEDIDENTITYALLREDUCE_H

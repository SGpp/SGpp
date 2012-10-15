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

		if (this->vecMode_ != OpenCL && this->vecMode_ != ArBB && this->vecMode_ != Hybrid_X86SIMD_OpenCL)
		{
			this->dataset_->transpose();
		}
		sg::parallel::PartitioningTool::getMPIPartitionSegment(numPatchedTrainingInstances_, &data_size, &data_offset, sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));

		this->level_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());
		this->index_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());

		m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
	}

	virtual void mult(base::DataVector &alpha, base::DataVector &result){
		sg::base::DataVector temp(this->numPatchedTrainingInstances_);
		sg::base::DataVector result_tmp(result.getSize());

		temp.setAll(0.0);
		result_tmp.setAll(0.0);
		result.setAll(0.0);

		size_t thread_count = 1;
		size_t thread_num = 0;
		this->myTimer_->start();
		#pragma omp parallel private(thread_num, thread_count)
		{
#ifdef _OPENMP
		thread_count = omp_get_num_threads();
		thread_num = omp_get_thread_num();
#endif
			size_t threadChunkStart, threadChunkEnd;
			sg::parallel::PartitioningTool::getPartitionSegment(
					data_offset, data_offset + data_size,
					thread_count, thread_num, &threadChunkStart, &threadChunkEnd,
					sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));
			MultType::mult(
					level_,
					index_,
					dataset_,
					alpha,
					temp,
					0,
					m_grid.getSize(),
					threadChunkStart,
					threadChunkEnd);
#pragma omp single nowait
			{
				double timeMult = this->myTimer_->stop();
				this->computeTimeMult_ += timeMult;
				this->completeTimeMult_ += timeMult;
				this->myTimer_->start();
			}

			// patch result -> set additional entries zero
			// only done for processes that need this part of the temp data
			if (data_offset+data_size > this->numTrainingInstances_ && this->numTrainingInstances_ != temp.getSize())
			{
				for (size_t i = 0; i < (temp.getSize()-this->numTrainingInstances_); i++)
				{
					temp.set(temp.getSize()-(i+1), 0.0f);
				}
			}

			sg::parallel::PartitioningTool::getPartitionSegment(
					0, m_grid.getSize(),
					thread_count, thread_num, &threadChunkStart, &threadChunkEnd, 1);
			MultTransType::multTranspose(
					level_,
					index_,
					dataset_,
					temp,
					result_tmp,
					threadChunkStart,
					threadChunkEnd,
					data_offset,
					data_offset + data_size);
		}
		this->computeTimeMultTrans_ += this->myTimer_->stop();
		MPI_Allreduce(result_tmp.getPointer(), result.getPointer(), result.getSize(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		this->completeTimeMultTrans_ += this->myTimer_->stop();

		if(sg::parallel::myGlobalMPIComm->getMyRank() == 0){
			std::cout << "*";
		}
	}

	virtual void generateb(base::DataVector &classes, base::DataVector &b){
		sg::base::DataVector myClasses(classes);
		base::DataVector b_tmp(b);
		b_tmp.setAll(0.0);
		b.setAll(0.0);
		// Apply padding
		if (this->numPatchedTrainingInstances_ != myClasses.getSize())
		{
			myClasses.resizeZero(this->numPatchedTrainingInstances_);
		}

		size_t thread_count = 1;
		size_t thread_num = 0;
		#pragma omp parallel private(thread_num, thread_count)
		{
#ifdef _OPENMP
		thread_count = omp_get_num_threads();
		thread_num = omp_get_thread_num();
#endif
			size_t threadChunkStart, threadChunkEnd;
			sg::parallel::PartitioningTool::getPartitionSegment(
					0, m_grid.getSize(),
					thread_count, thread_num, &threadChunkStart, &threadChunkEnd, 1);
			MultTransType::multTranspose(
						level_,
						index_,
						dataset_,
						myClasses,
						b_tmp,
						threadChunkStart,
						threadChunkEnd,
						data_offset,
						data_offset + data_size);
		}
		MPI_Allreduce(b_tmp.getPointer(), b.getPointer(), b_tmp.getSize(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}

	virtual void rebuildLevelAndIndex(){
		delete this->level_;
		delete this->index_;

		this->level_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());
		this->index_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());

		m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
	}

};

}
}

#endif // DMSYSTEMMATRIXVECTORIZEDIDENTITYALLREDUCE_H

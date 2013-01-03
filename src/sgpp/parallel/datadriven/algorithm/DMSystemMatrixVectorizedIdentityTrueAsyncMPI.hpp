/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Roman Karlstetter (karlstetter@mytum.de)
#ifndef DMSYSTEMMATRIXVECTORIZEDIDENTITYTRUEASYNCMPI_HPP
#define DMSYSTEMMATRIXVECTORIZEDIDENTITYTRUEASYNCMPI_HPP


#include "datadriven/algorithm/DMSystemMatrixBase.hpp"
#include "parallel/tools/MPI/SGppMPITools.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/exception/operation_exception.hpp"
#include "base/grid/Grid.hpp"
#include "parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "parallel/tools/PartitioningTool.hpp"
#include "parallel/tools/TypesParallel.hpp"

#include <strstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sg
{
namespace parallel
{

template<typename MultType, typename MultTransType>
class DMSystemMatrixVectorizedIdentityTrueAsyncMPI : public sg::datadriven::DMSystemMatrixBase
{
private:
	/// vectorization mode
	VectorizationType vecMode_;
	/// Number of orignal training instances
	size_t numTrainingInstances_;
	/// Number of patched and used training instances
	size_t numPatchedTrainingInstances_;

	sg::base::DataMatrix* level_;
	/// Member to store the sparse grid's indices for better vectorization
	sg::base::DataMatrix* index_;

	/// reference to grid. needed to get new grid size after it changes
	sg::base::Grid& m_grid;

	/// how to distribute storage array across processes
	int* _mpi_grid_sizes;
	int* _mpi_grid_offsets;

	/// how to distribute dataset across processes
	int* _mpi_data_sizes;
	int* _mpi_data_offsets;

public:
	/**
	 * Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param trainData reference to sg::base::DataMatrix that contains the training data
	 * @param lambda the lambda, the regression parameter
	 * @param vecMode vectorization mode
	 */
	DMSystemMatrixVectorizedIdentityTrueAsyncMPI(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
		: DMSystemMatrixBase(trainData, lambda), vecMode_(vecMode), numTrainingInstances_(0), numPatchedTrainingInstances_(0), m_grid(SparseGrid)
	{
		// handle unsupported vector extensions
		if (this->vecMode_ != X86SIMD &&
				this->vecMode_ != MIC &&
				this->vecMode_ != Hybrid_X86SIMD_MIC &&
				this->vecMode_ != OpenCL &&
				this->vecMode_ != ArBB &&
				this->vecMode_ != Hybrid_X86SIMD_OpenCL)
		{
			throw new sg::base::operation_exception("DMSystemMatrixSPVectorizedIdentity : un-supported vector extension!");
		}

		this->dataset_ = new sg::base::DataMatrix(trainData);
		this->numTrainingInstances_ = this->dataset_->getNrows();
		this->numPatchedTrainingInstances_ = sg::parallel::DMVectorizationPaddingAssistant::padDataset(*(this->dataset_), vecMode_);

		if (this->vecMode_ != OpenCL && this->vecMode_ != ArBB && this->vecMode_ != Hybrid_X86SIMD_OpenCL)
		{
			this->dataset_->transpose();
		}

		int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();

		// arrays for distribution settings
		_mpi_data_sizes = new int[mpi_size];
		_mpi_data_offsets = new int[mpi_size];

		sg::parallel::PartitioningTool::calcDistribution(this->numPatchedTrainingInstances_, mpi_size, _mpi_data_sizes, _mpi_data_offsets, sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));

		// arrays for distribution settings
		_mpi_grid_sizes = new int[mpi_size];
		_mpi_grid_offsets = new int[mpi_size];
		sg::parallel::PartitioningTool::calcDistribution(m_grid.getSize(), mpi_size, _mpi_grid_sizes, _mpi_grid_offsets, 1);

		if(sg::parallel::myGlobalMPIComm->getMyRank() == 0){
			std::cout << "Data: " << mpi_size << " chunks of max size " << _mpi_data_sizes[0] << " (blocksize:" << sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_) << "), total: " << this->numPatchedTrainingInstances_ << std::endl;
			std::cout << "Grid: " << mpi_size << " chunks of max size " << _mpi_grid_sizes[0] << " (blocksize: 1), total: " << m_grid.getSize() << std::endl;
		}
		this->level_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());
		this->index_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());

		m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

		// mult: distribute calculations over dataset
		// multTranspose: distribute calculations over grid
	}

	/**
	 * Std-Destructor
	 */
	virtual ~DMSystemMatrixVectorizedIdentityTrueAsyncMPI(){
		delete this->dataset_;

		delete this->level_;
		delete this->index_;

		delete[] this->_mpi_grid_sizes;
		delete[] this->_mpi_grid_offsets;
		delete[] this->_mpi_data_sizes;
		delete[] this->_mpi_data_offsets;
	}

	virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result){
		sg::base::DataVector temp(this->numPatchedTrainingInstances_);
		//sg::base::DataVector temp_reference(this->numPatchedTrainingInstances_);
		sg::base::DataVector result_reference(result.getSize());

		result.setAll(0.0);
		temp.setAll(0.0);
		//temp_reference.setAll(0.0);
		result_reference.setAll(0.0);
		double* ptrResult = result.getPointer();
		double* ptrTemp = temp.getPointer();

		int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
		int mpi_myrank = sg::parallel::myGlobalMPIComm->getMyRank();

		//MultType::mult(level_, index_, dataset_, alpha, temp_reference, 0, alpha.getSize(), 0, this->numPatchedTrainingInstances_);
		for (size_t i = this->numTrainingInstances_; i < this->numPatchedTrainingInstances_; i++)
		{
			//temp_reference.set(i, 0.0f);
		}

		MPI_Request dataSendReqs[mpi_size];
		MPI_Request dataRecvReqs[mpi_size]; //allocating a little more than necessary, otherwise complicated index computations needed
		int tagsData[mpi_size];
		for (int i = 0; i<mpi_size; i++){
			tagsData[i] = _mpi_data_offsets[i]*2 + 2;
		}
		for(int rank = 0; rank <mpi_size; rank++){
			if(rank == mpi_myrank){
				dataRecvReqs[rank] = MPI_REQUEST_NULL;
				continue;
			}
			MPI_Irecv(&ptrTemp[_mpi_data_offsets[rank]], _mpi_data_sizes[rank], MPI_DOUBLE, rank, tagsData[rank], MPI_COMM_WORLD, &dataRecvReqs[rank]);
		}

		MPI_Request gridSendReqs[mpi_size];
		MPI_Request gridRecvReqs[mpi_size]; //allocating a little more than necessary, otherwise complicated index computations needed
		int tagsGrid[mpi_size];
		for (int i = 0; i<mpi_size; i++){
			tagsGrid[i] = _mpi_grid_offsets[i]*2 + 3;
		}
		for(int rank = 0; rank <mpi_size; rank++){
			if(rank == mpi_myrank){
				gridRecvReqs[rank] = MPI_REQUEST_NULL;
				continue;
			}
			MPI_Irecv(&ptrResult[_mpi_grid_offsets[rank]], _mpi_grid_sizes[rank], MPI_DOUBLE, rank, tagsGrid[rank], MPI_COMM_WORLD, &gridRecvReqs[rank]);
		}
//		MPI_Request dataSendReqs[_mpi_data_sizes_global[mpi_myrank] * mpi_size];
//		MPI_Request gridSendReqs[_mpi_grid_sizes_global[mpi_myrank] * mpi_size];

		this->myTimer_->start();

		size_t dataProcessChunkStart = _mpi_data_offsets[mpi_myrank];
		size_t dataProcessChunkEnd = dataProcessChunkStart + _mpi_data_sizes[mpi_myrank];
		size_t gridProcessChunkStart = _mpi_grid_offsets[mpi_myrank];
		size_t gridProcessChunkEnd = gridProcessChunkStart + _mpi_grid_sizes[mpi_myrank];

	#pragma omp parallel
		{
			size_t threadStartData, threadEndData;

			sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
					dataProcessChunkStart, dataProcessChunkEnd,
					&threadStartData, &threadEndData, 1);

			MultType::mult(level_, index_, dataset_, alpha, temp, 0, alpha.getSize(), threadStartData, threadEndData);
				// patch result -> set additional entries zero
				// only done for processes that need this part of the temp data for multTrans
			for (size_t i = std::max<size_t>(this->numTrainingInstances_, threadStartData); i < threadEndData; i++)
			{
				temp.set(i, 0.0f);
			}
//				for(int i = start; i<end;i++){
//					if(temp[i] != temp_reference[i]){
//						std::cout << "results do not match: " << i << ": " <<temp[i] << "!=" << temp_reference[i] << "thread/proc" << thread_num <<"/" << mpi_myrank << std::endl;
//						throw new sg::base::operation_exception("reuslts do not match!");
//					}
//				}


#pragma omp barrier // make sure that all threads finished their part
#pragma omp single nowait
			{
				sg::parallel::myGlobalMPIComm->IsendToAll(&ptrTemp[dataProcessChunkStart], _mpi_data_sizes[mpi_myrank],
														  dataProcessChunkStart*2 + 2, &dataSendReqs[mpi_myrank]);
			}

			size_t threadStartGrid, threadEndGrid;

			sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
				gridProcessChunkStart, gridProcessChunkEnd,
				&threadStartGrid, &threadEndGrid, 1);

			MultTransType::multTranspose(
					level_,
					index_,
					dataset_,
					temp,
					result,
					threadStartGrid,
					threadEndGrid,
					dataProcessChunkStart,
					dataProcessChunkEnd);
			while (true) {
				int idx;
				MPI_Waitany(mpi_size, dataRecvReqs, &idx, MPI_STATUS_IGNORE); // MPI_Waitany is thread safe
				if(idx == MPI_UNDEFINED) {
					// no more active request, everything is done
					break;
				}
				size_t dataChunkStart = _mpi_data_offsets[idx];
				size_t dataChunkEnd = dataChunkStart + _mpi_data_sizes[idx];

				MultTransType::multTranspose(
							level_,
							index_,
							dataset_,
							temp,
							result,
							threadStartGrid,
							threadEndGrid,
							dataChunkStart,
							dataChunkEnd);
			}
		#pragma omp barrier // wait for all threads to finish

		}
		sg::parallel::myGlobalMPIComm->IsendToAll(&ptrResult[gridProcessChunkStart], _mpi_data_sizes[mpi_myrank],
												  gridProcessChunkStart*2 + 3, &gridSendReqs[mpi_myrank]);

		this->computeTimeMultTrans_ += this->myTimer_->stop();

		MPI_Waitall(mpi_size, gridRecvReqs, MPI_STATUSES_IGNORE);
		this->completeTimeMultTrans_ += this->myTimer_->stop();

		result.axpy(static_cast<double>(this->numTrainingInstances_)*this->lambda_, alpha);
		if (mpi_myrank == 0) std::cout << "*";
	}

	virtual void generateb(sg::base::DataVector& classes, sg::base::DataVector& b){
		int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
		int mpi_myrank = sg::parallel::myGlobalMPIComm->getMyRank();

		double* ptrB = b.getPointer();
		b.setAll(0.0);

		sg::base::DataVector myClasses(classes);
		// Apply padding
		if (this->numPatchedTrainingInstances_ != myClasses.getSize())
		{
			myClasses.resizeZero(this->numPatchedTrainingInstances_);
		}


	#pragma omp parallel
		{
			size_t threadChunkStart;
			size_t threadChunkEnd;
			sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
					_mpi_grid_offsets[mpi_myrank], _mpi_grid_offsets[mpi_myrank] + _mpi_grid_sizes[mpi_myrank],
					&threadChunkStart, &threadChunkEnd, 1);
			MultTransType::multTranspose(level_, index_, dataset_, myClasses, b, threadChunkStart, threadChunkEnd, 0, this->numPatchedTrainingInstances_);
	#pragma omp barrier
		}
		sg::parallel::myGlobalMPIComm->dataVectorAllToAll(b, _mpi_grid_offsets, _mpi_grid_sizes);
	}

	virtual void rebuildLevelAndIndex(){
		delete this->level_;
		delete this->index_;

		this->level_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());
		this->index_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());

		m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

		int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();

		sg::parallel::PartitioningTool::calcDistribution(m_grid.getSize(), mpi_size, _mpi_grid_sizes, _mpi_grid_offsets, 1);
	}


};

}
}

#endif // DMSYSTEMMATRIXVECTORIZEDIDENTITYTRUEASYNCMPI_HPP

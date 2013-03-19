/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef DMSYSTEMMATRIXVECTORIZEDIDENTITYASYNCMPI_HPP
#define DMSYSTEMMATRIXVECTORIZEDIDENTITYASYNCMPI_HPP

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

/**
 * Class that implements the virtual class sg::base::OperationMatrix for the
 * application of classification for the Systemmatrix
 *
 * The Identity matrix is used as regularization operator.
 *
 * For the Operation B's mult and mutlTransposed functions
 * vectorized formulations are used.
 */
template<typename MultType, typename MultTransType>
class DMSystemMatrixVectorizedIdentityAsyncMPI : public sg::datadriven::DMSystemMatrixBase
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


public:
	/**
	 * Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param trainData reference to sg::base::DataMatrix that contains the training data
	 * @param lambda the lambda, the regression parameter
	 * @param vecMode vectorization mode
	 */
	DMSystemMatrixVectorizedIdentityAsyncMPI(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
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


		/* calculate distribution of data */
		calculateChunkCountData();
		_mpi_data_sizes = new int[_chunkCountData];
		_mpi_data_offsets = new int[_chunkCountData];
		_mpi_data_sizes_global = new int[mpi_size];
		_mpi_data_offsets_global = new int[mpi_size];
		sg::parallel::PartitioningTool::calcDistribution(this->numPatchedTrainingInstances_, _chunkCountData, _mpi_data_sizes, _mpi_data_offsets, sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));
		sg::parallel::PartitioningTool::calcDistribution(_chunkCountData, mpi_size, _mpi_data_sizes_global, _mpi_data_offsets_global, 1);
		if(sg::parallel::myGlobalMPIComm->getMyRank() == 0){
			std::cout << "Max size per chunk Data: " << _mpi_data_sizes[0] << std::endl;
		}

		/* calculate distribution of grid */
		calculateChunkCountGrid();
		_mpi_grid_sizes = new int[_chunkCountGrid];
		_mpi_grid_offsets = new int[_chunkCountGrid];
		_mpi_grid_sizes_global = new int[mpi_size];
		_mpi_grid_offsets_global = new int[mpi_size];
		sg::parallel::PartitioningTool::calcDistribution(m_grid.getSize(), _chunkCountGrid, _mpi_grid_sizes, _mpi_grid_offsets, 1);
		sg::parallel::PartitioningTool::calcDistribution(_chunkCountGrid, mpi_size, _mpi_grid_sizes_global, _mpi_grid_offsets_global, 1);


		this->level_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());
		this->index_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());
		m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

		// mult: distribute calculations over dataset
		// multTranspose: distribute calculations over grid
	}

	/**
	 * Destructor
	 */
	virtual ~DMSystemMatrixVectorizedIdentityAsyncMPI(){
		delete this->dataset_;

		delete this->level_;
		delete this->index_;

		delete[] this->_mpi_grid_sizes;
		delete[] this->_mpi_grid_offsets;
		delete[] this->_mpi_data_sizes;
		delete[] this->_mpi_data_offsets;

		delete[] this->_mpi_grid_sizes_global;
		delete[] this->_mpi_grid_offsets_global;
		delete[] this->_mpi_data_sizes_global;
		delete[] this->_mpi_data_offsets_global;
	}

	virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result){

		sg::base::DataVector temp(this->numPatchedTrainingInstances_);
		result.setAll(0.0);
		temp.setAll(0.0);
		double* ptrResult = result.getPointer();
		double* ptrTemp = temp.getPointer();

		int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
		int mpi_myrank = sg::parallel::myGlobalMPIComm->getMyRank();

		/* setup MPI_Requests, tags and post receives for data */
		MPI_Request dataRecvReqs[_chunkCountData]; //allocating a little more than necessary, otherwise complicated index computations needed
		int tagsData[_chunkCountData];
		for (int i = 0; i<_chunkCountData; i++){
			tagsData[i] = _mpi_data_offsets[i]*2 + 2;
		}
		sg::parallel::myGlobalMPIComm->IrecvFromAll(ptrTemp, _mpi_data_sizes_global, _mpi_data_offsets_global, _mpi_data_sizes, _mpi_data_offsets, tagsData, _chunkCountData, dataRecvReqs);

		/* setup MPI_Requests, tags and post receives for grid */
		MPI_Request gridRecvReqs[_chunkCountGrid]; //allocating a little more than necessary, otherwise complicated index computations needed
		int tagsGrid[_chunkCountGrid];
		for (int i = 0; i<_chunkCountGrid; i++){
			tagsGrid[i] = _mpi_grid_offsets[i]*2 + 3;
		}
		sg::parallel::myGlobalMPIComm->IrecvFromAll(ptrResult, _mpi_grid_sizes_global, _mpi_grid_offsets_global, _mpi_grid_sizes, _mpi_grid_offsets, tagsGrid, _chunkCountGrid, gridRecvReqs);
		MPI_Request dataSendReqs[_mpi_data_sizes_global[mpi_myrank] * mpi_size];
		MPI_Request gridSendReqs[_mpi_grid_sizes_global[mpi_myrank] * mpi_size];

		this->myTimer_->start();
		#pragma omp parallel
		{
			size_t myDataChunkStart = _mpi_data_offsets_global[mpi_myrank];
			size_t myDataChunkEnd = myDataChunkStart + _mpi_data_sizes_global[mpi_myrank];
			size_t threadStart, threadEnd;
			sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
					myDataChunkStart, myDataChunkEnd,
					&threadStart, &threadEnd, 1);
			for(size_t chunkIndex = threadStart; chunkIndex < threadEnd; chunkIndex++){
				size_t start = _mpi_data_offsets[chunkIndex];
				size_t end =  start + _mpi_data_sizes[chunkIndex];
				MultType::mult(level_, index_, dataset_, alpha, temp, 0, alpha.getSize(), start, end);
				// patch result -> set additional entries zero
				// only done for processes that need this part of the temp data for multTrans
				for (size_t i = std::max<size_t>(this->numTrainingInstances_, start); i < end; i++)
				{
					temp.set(i, 0.0f);
				}
				sg::parallel::myGlobalMPIComm->IsendToAll(&ptrTemp[start], _mpi_data_sizes[chunkIndex], tagsData[chunkIndex], &dataSendReqs[(chunkIndex - myDataChunkStart)*mpi_size]);
			}

#pragma omp single
			{
				double computationTime = this->myTimer_->stop();
				this->computeTimeMult_ += computationTime;
				if(MPI_Waitall(_chunkCountData, dataRecvReqs, MPI_STATUSES_IGNORE) != MPI_SUCCESS){
					std::cout << "errors in communication" << std::endl;
				}
				// we don't really need to wait for the sends to
				// finish as we don't need (in particular not modify) temp
				// advantage: it's a lot faster like this
//				if(MPI_Waitall(_mpi_data_sizes_global[mpi_myrank] * mpi_size, dataSendReqs, MPI_STATUSES_IGNORE) != MPI_SUCCESS){
//					std::cout << "errors in communication (send)" << std::endl;
//				}
				double completeTime = this->myTimer_->stop();
				this->completeTimeMult_ += completeTime;

//				double communicationTime = completeTime - computationTime;
//				double maxComputationTime, minComputationTime;
//				double maxCompleteTime, minCompleteTime;
//				double maxCommunicationTime, minCommunicationTime;

//				MPI_Reduce(&computationTime, &maxComputationTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
//				MPI_Reduce(&computationTime, &minComputationTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
//				MPI_Reduce(&communicationTime, &maxCommunicationTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
//				MPI_Reduce(&communicationTime, &minCommunicationTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
//				MPI_Reduce(&completeTime, &maxCompleteTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
//				MPI_Reduce(&completeTime, &minCompleteTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
//				if (sg::parallel::myGlobalMPIComm->getMyRank() == 0 && false) {
//					std::cout << "size: " << _mpi_data_sizes[0]*sizeof(double) << "(sizeof (double): "<< sizeof(double) <<")" << std::endl;
//					std::cout << "computation     time min - max: " << minComputationTime << " - " << maxComputationTime << " (difference: " << (maxComputationTime - minComputationTime) << ") " << std::endl;
//					std::cout << "complete        time min - max: " << minCompleteTime << " - " << maxCompleteTime << " (difference: " << (maxCompleteTime - minCompleteTime) << ") " << std::endl;
//					std::cout << "communication   time min - max: " << minCommunicationTime << " - " << maxCommunicationTime << " (difference: " << (maxCommunicationTime - minCommunicationTime) << ") " << std::endl;
//					std::cout << std::endl;
//				}
//				/// @todo min/max mpireduce to see where the time is spent

				this->myTimer_->start();

			} //implicit OpenMP barrier here

			size_t myGridChunkStart = _mpi_grid_offsets_global[mpi_myrank];
			size_t myGridChunkEnd = myGridChunkStart + _mpi_grid_sizes_global[mpi_myrank];
			sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
					myGridChunkStart, myGridChunkEnd, &threadStart, &threadEnd, 1);

			for(size_t thread_chunk = threadStart; thread_chunk<threadEnd; thread_chunk++){
				size_t start = _mpi_grid_offsets[thread_chunk];
				size_t end =  start + _mpi_grid_sizes[thread_chunk];

				MultTransType::multTranspose(level_, index_, dataset_, temp, result, start, end, 0, this->numPatchedTrainingInstances_);
				sg::parallel::myGlobalMPIComm->IsendToAll(&ptrResult[start], _mpi_grid_sizes[thread_chunk], tagsGrid[thread_chunk], &gridSendReqs[(thread_chunk - myGridChunkStart)*mpi_size]);
			}
		}
		this->computeTimeMultTrans_ += this->myTimer_->stop();

		MPI_Status gridRecsStats[_chunkCountGrid]; //allocating a little more than necessary, otherwise complicated index computations needed
		if(MPI_Waitall(_chunkCountGrid, gridRecvReqs, gridRecsStats) != MPI_SUCCESS) {
			std::cout << "Communication error (waitall gridrecvreqs)" << std::endl;
			throw new sg::base::operation_exception("Communication error!");
		}
		if(MPI_Waitall(_chunkCountGrid, gridSendReqs, MPI_STATUSES_IGNORE) != MPI_SUCCESS) {
			std::cout << "Communication error (waitall gridsendreqs)" << std::endl;
			throw new sg::base::operation_exception("Communication error!");
		}
		this->completeTimeMultTrans_ += this->myTimer_->stop();

		result.axpy(static_cast<double>(this->numTrainingInstances_)*this->lambda_, alpha);
		if (mpi_myrank == 0) std::cout << "*";
	} //end mult

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

		MPI_Request gridRecvReqs[_chunkCountGrid]; //allocating a little more than necessary, otherwise complicated index computations needed
		int tags[_chunkCountGrid];
		for(int i = 0; i<_chunkCountGrid; i++){
			tags[i] = _mpi_grid_offsets[i]+1;
		}
		sg::parallel::myGlobalMPIComm->IrecvFromAll(ptrB, _mpi_grid_sizes_global, _mpi_grid_offsets_global, _mpi_grid_sizes, _mpi_grid_offsets, tags, _chunkCountGrid, gridRecvReqs);
		MPI_Request gridSendReqs[mpi_size * _mpi_grid_sizes_global[mpi_myrank]];
#pragma omp parallel
		{
			size_t myGridChunkStart = _mpi_grid_offsets_global[mpi_myrank];
			size_t myGridChunkEnd = myGridChunkStart + _mpi_grid_sizes_global[mpi_myrank];
			size_t threadStart, threadEnd;
			sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
					myGridChunkStart, myGridChunkEnd, &threadStart, &threadEnd, 1);

			for(size_t thread_chunk = threadStart; thread_chunk<threadEnd; thread_chunk++){
				size_t start = _mpi_grid_offsets[thread_chunk];
				size_t end =  start + _mpi_grid_sizes[thread_chunk];

				MultTransType::multTranspose(level_, index_, dataset_, myClasses, b, start, end, 0, this->numPatchedTrainingInstances_);
				sg::parallel::myGlobalMPIComm->IsendToAll(&ptrB[start], _mpi_grid_sizes[thread_chunk], start+1, &gridSendReqs[(thread_chunk - myGridChunkStart)*mpi_size]);
			}
		}
		MPI_Status stats[_chunkCountGrid];
		if(MPI_Waitall(_chunkCountGrid, gridRecvReqs, stats) != MPI_SUCCESS) {
			std::cout << "communication error (recvReqs) in generateB" << std::endl;
			throw new sg::base::operation_exception("COmmunication Error");
		}
		if(MPI_Waitall(mpi_size * _mpi_grid_sizes_global[mpi_myrank], gridSendReqs, stats) != MPI_SUCCESS) {
			std::cout << "communication error (sendReqs) in generateB" << std::endl;
			throw new sg::base::operation_exception("COmmunication Error");
		}
	}

	virtual void rebuildLevelAndIndex(){
		delete this->level_;
		delete this->index_;
		this->level_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());
		this->index_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());

		m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));


		delete[] _mpi_grid_sizes;
		delete[] _mpi_grid_offsets;

		calculateChunkCountGrid();
		_mpi_grid_sizes = new int[_chunkCountGrid];
		_mpi_grid_offsets = new int[_chunkCountGrid];
		sg::parallel::PartitioningTool::calcDistribution(m_grid.getSize(), _chunkCountGrid, _mpi_grid_sizes, _mpi_grid_offsets, 1);
		sg::parallel::PartitioningTool::calcDistribution(_chunkCountGrid, sg::parallel::myGlobalMPIComm->getNumRanks(), _mpi_grid_sizes_global, _mpi_grid_offsets_global, 1);
	}

private:
	/// reference to grid. needed to get new grid size after it changes
	sg::base::Grid& m_grid;

	/// how to distribute storage array across processes
	int * _mpi_grid_sizes;
	int * _mpi_grid_offsets;

	/// how to distribute dataset across processes
	int * _mpi_data_sizes;
	int * _mpi_data_offsets;

	/// which chunks belong to which process
	int * _mpi_data_sizes_global;
	int * _mpi_data_offsets_global;

	/// which chunks belong to which process
	int * _mpi_grid_sizes_global;
	int * _mpi_grid_offsets_global;

	/// into how many chunks should data and grid be partitioned
	size_t _chunkCountData;
	size_t _chunkCountGrid;

#define APPROXCHUNKSIZEGRID_ASYNC 20 // approximately how many blocks should be computed for the grid before sending
#define APPROXCHUNKSIZEDATA_ASYNC 150 // approximately how many blocks should be computed for the data before sending

	void calculateChunkCountGrid() {
		_chunkCountGrid = m_grid.getSize()/APPROXCHUNKSIZEGRID_ASYNC;
		_chunkCountGrid = std::max<size_t>(_chunkCountGrid, getMinChunkCount());
		int chunkSize =
				static_cast<int>(
					std::ceil(
						static_cast<double>(m_grid.getSize())/static_cast<double>(_chunkCountGrid)
						)
					);
		if(sg::parallel::myGlobalMPIComm->getMyRank() == 0){
			std::cout << "Number of chunks Grid: " << _chunkCountGrid << std::endl;
			std::cout << "Max size per chunk Grid: " << chunkSize << std::endl;
		}
	}

	void calculateChunkCountData() {
		size_t blockCount = this->numPatchedTrainingInstances_/sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_); // this process has (in total) blockCountData blocks of Data to process
		_chunkCountData = blockCount/APPROXCHUNKSIZEDATA_ASYNC;
		_chunkCountData = std::max<size_t>(_chunkCountData, getMinChunkCount());
		if(sg::parallel::myGlobalMPIComm->getMyRank() == 0){
			std::cout << "Number of chunks Data: " << _chunkCountData << std::endl;
		}
	}

	size_t getMinChunkCount(){
		int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
		size_t thread_count = 1;
#ifdef _OPENMP
#pragma omp parallel
		{
			thread_count = omp_get_num_threads();
		}
#endif

		// every process should have at least thread_count chunks,
		// such that every thread has at least one chunk
		return mpi_size * thread_count;

	}
};

}
}

#endif /* DMSYSTEMMATRIXVECTORIZEDIDENTITYASYNCMPI_HPP */

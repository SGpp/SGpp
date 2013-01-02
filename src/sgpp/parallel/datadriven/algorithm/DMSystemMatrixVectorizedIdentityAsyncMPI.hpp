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

#define APPROXCHUNKSIZEGRID_ASYNC 100 // approximately how many blocks should be computed for the grid before sending
#define APPROXCHUNKSIZEDATA_ASYNC 100 // approximately how many blocks should be computed for the data before sending

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
		size_t thread_count = 1;
#ifdef _OPENMP
#pragma omp parallel
		{
			thread_count = omp_get_num_threads();
		}
#endif

		size_t minChunkCount = mpi_size * thread_count; // every process should have at least thread_count chunks

		size_t blockCountData = this->numPatchedTrainingInstances_/sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_); // this process has (in total) blockCountData blocks of Data to process
		_chunkCountData = blockCountData/APPROXCHUNKSIZEDATA_ASYNC;
		_chunkCountData = std::max<size_t>(_chunkCountData, minChunkCount);

		// arrays for distribution settings
		_mpi_data_sizes = new int[_chunkCountData];
		_mpi_data_offsets = new int[_chunkCountData];
		_mpi_data_sizes_global = new int[mpi_size];
		_mpi_data_offsets_global = new int[mpi_size];

		sg::parallel::PartitioningTool::calcDistribution(this->numPatchedTrainingInstances_, _chunkCountData, _mpi_data_sizes, _mpi_data_offsets, sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));
		sg::parallel::PartitioningTool::calcDistribution(_chunkCountData, mpi_size, _mpi_data_sizes_global, _mpi_data_offsets_global, 1);

		_chunkCountGrid = m_grid.getSize()/APPROXCHUNKSIZEGRID_ASYNC;
		_chunkCountGrid = std::max<size_t>(_chunkCountGrid, minChunkCount); // every process should have at least thread_count chunk

		// arrays for distribution settings
		_mpi_grid_sizes = new int[_chunkCountGrid];
		_mpi_grid_offsets = new int[_chunkCountGrid];
		_mpi_grid_sizes_global = new int[mpi_size];
		_mpi_grid_offsets_global = new int[mpi_size];

		sg::parallel::PartitioningTool::calcDistribution(m_grid.getSize(), _chunkCountGrid, _mpi_grid_sizes, _mpi_grid_offsets, 1);
		sg::parallel::PartitioningTool::calcDistribution(_chunkCountGrid, mpi_size, _mpi_grid_sizes_global, _mpi_grid_offsets_global, 1);

		if(sg::parallel::myGlobalMPIComm->getMyRank() == 0){
			std::cout << "Data: " << _chunkCountData << " chunks of max size " << _mpi_data_sizes[0] << " (blocksize:" << sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_) << "), total: " << this->numPatchedTrainingInstances_ << std::endl;
			std::cout << "Grid: " << _chunkCountGrid << " chunks of max size " << _mpi_grid_sizes[0] << " (blocksize: 1), total: " << m_grid.getSize() << std::endl;
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
		MPI_Request dataRecvReqs[_chunkCountData]; //allocating a little more than necessary, otherwise complicated index computations needed
		int tagsData[_chunkCountData];
		for (int i = 0; i<_chunkCountData; i++){
			tagsData[i] = _mpi_data_offsets[i]*2 + 2;
		}
		sg::parallel::myGlobalMPIComm->IrecvFromAll(ptrTemp, _mpi_data_sizes_global, _mpi_data_offsets_global, _mpi_data_sizes, _mpi_data_offsets, tagsData, _chunkCountData, dataRecvReqs);

		MPI_Request gridRecvReqs[_chunkCountGrid]; //allocating a little more than necessary, otherwise complicated index computations needed
		int tagsGrid[_chunkCountGrid];
		for (int i = 0; i<_chunkCountGrid; i++){
			tagsGrid[i] = _mpi_grid_offsets[i]*2 + 3;
		}
		sg::parallel::myGlobalMPIComm->IrecvFromAll(ptrResult, _mpi_grid_sizes_global, _mpi_grid_offsets_global, _mpi_grid_sizes, _mpi_grid_offsets, tagsGrid, _chunkCountGrid, gridRecvReqs);

		MPI_Request dataSendReqs[_mpi_data_sizes_global[mpi_myrank] * mpi_size];
		MPI_Request gridSendReqs[_mpi_grid_sizes_global[mpi_myrank] * mpi_size];

		this->myTimer_->start();
		int thread_count = 1;
		int thread_num = 0;
		#pragma omp parallel private (thread_num, thread_count)
		{
	#ifdef _OPENMP
			thread_count = omp_get_num_threads();
			thread_num = omp_get_thread_num();
	#endif
			size_t myDataChunkStart = _mpi_data_offsets_global[mpi_myrank];
			size_t myDataChunkEnd = myDataChunkStart + _mpi_data_sizes_global[mpi_myrank];
			size_t threadStart, threadEnd;
			sg::parallel::PartitioningTool::getPartitionSegment(
					myDataChunkStart, myDataChunkEnd,
					thread_count, thread_num, &threadStart, &threadEnd, 1);
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
//				for(int i = start; i<end;i++){
//					if(temp[i] != temp_reference[i]){
//						std::cout << "results do not match: " << i << ": " <<temp[i] << "!=" << temp_reference[i] << "thread/proc" << thread_num <<"/" << mpi_myrank << std::endl;
//						throw new sg::base::operation_exception("reuslts do not match!");
//					}
//				}
				sg::parallel::myGlobalMPIComm->IsendToAll(&ptrTemp[start], _mpi_data_sizes[chunkIndex], start*2 + 2, &dataSendReqs[(chunkIndex - myDataChunkStart)*mpi_size]);
			}
//#pragma omp barrier
//			for(int chunkIndex = myDataChunkStart; chunkIndex<myDataChunkEnd;chunkIndex++){
//				for(int i = _mpi_data_offsets[chunkIndex]; i<_mpi_data_offsets[chunkIndex] + _mpi_data_sizes[chunkIndex];i++){
//					if(temp[i] != temp_reference[i]){
//						std::cout << "results do not match (after barrier): " << i << ": " <<temp[i] << "!=" << temp_reference[i] << "thread/proc" << thread_num <<"/" << mpi_myrank << std::endl;
//						throw new sg::base::operation_exception("reuslts do not match!");
//					}
//				}
//			}
//			for(size_t chunkIndex = threadStart; chunkIndex < threadEnd; chunkIndex++){
//				size_t start = _mpi_data_offsets[chunkIndex];
//				sg::parallel::myGlobalMPIComm->IsendToAll(&ptrTemp[start], _mpi_data_sizes[chunkIndex], start*2 + 2, &dataSendReqs[(chunkIndex - myDataChunkStart)*mpi_size]);
//			}
//			for(int chunkIndex = myDataChunkStart; chunkIndex<myDataChunkEnd;chunkIndex++){
//				for(int i = _mpi_data_offsets[chunkIndex]; i<_mpi_data_offsets[chunkIndex] + _mpi_data_sizes[chunkIndex];i++){
//					if(temp[i] != temp_reference[i]){
//						std::cout << "results do not match (after send): " << i << ": " <<temp[i] << "!=" << temp_reference[i] << "thread/proc" << thread_num <<"/" << mpi_myrank << std::endl;
//						throw new sg::base::operation_exception("reuslts do not match!");
//					}
//				}
//			}

//			//for(size_t thread_chunk = myDataChunkStart + thread_num; thread_chunk < myDataChunkEnd; thread_chunk+=thread_count){
//			for(size_t thread_chunk = myDataChunkStart; thread_chunk < myDataChunkEnd; thread_chunk++){
//				if((thread_chunk - thread_num) % thread_count != 0){
//					continue;
//				}
//				size_t start = _mpi_data_offsets[thread_chunk];
//				size_t end =  start + _mpi_data_sizes[thread_chunk];

//				if (start % MultType::getChunkDataPoints() != 0 || end % MultType::getChunkDataPoints() != 0)
//				{
//					std::cout << "start%CHUNKDATAPOINTS_X86: " << start%sg::parallel::X86SimdLinearMult::getChunkDataPoints() << "; end%CHUNKDATAPOINTS_X86: " << end%sg::parallel::X86SimdLinearMult::getChunkDataPoints() << std::endl;
//					throw sg::base::operation_exception("processed vector segment must fit to CHUNKDATAPOINTS_X86!");
//				}
////				if(sg::parallel::myGlobalMPIComm->getMyRank() == 0){
////						std::cout << "threadchunk: " << thread_chunk <<  std::endl;
////				}


//				MultType::mult(level_, index_, dataset_, alpha, temp, 0, alpha.getSize(), start, end);


//				// patch result -> set additional entries zero
//				// only done for processes that need this part of the temp data for multTrans
//				for (size_t i = std::max<size_t>(this->numTrainingInstances_, start); i < end; i++)
//				{
//					temp.set(i, 0.0f);
//				}
//				for(int i = start; i<end;i++){
//					if(temp[i] != temp_reference[i]){
//						std::cout << "results do not match: " << i << ": " <<temp[i] << "!=" << temp_reference[i] << "thread/proc" << thread_num <<"/" << mpi_myrank << std::endl;
//						throw new sg::base::operation_exception("reuslts do not match!");
//					}
//				}
//				sg::parallel::myGlobalMPIComm->IsendToAll(&ptrTemp[start], _mpi_data_sizes[thread_chunk], start*2 + 2, &dataSendReqs[(thread_chunk - myDataChunkStart)*mpi_size]);
//			}


#pragma omp barrier //don't know why, but this line makes it fast...
#pragma omp single
			{
				double computationTime = this->myTimer_->stop();
				this->computeTimeMult_ += computationTime;
				MPI_Status dataRecvStats[_chunkCountData]; //allocating a little more than necessary, otherwise complicated index computations needed
				if(MPI_Waitall(_chunkCountData, dataRecvReqs, dataRecvStats) != MPI_SUCCESS){
					std::cout << "errors in communication" << std::endl;
				}
//				if(MPI_Waitall(_mpi_data_sizes_global[mpi_myrank] * mpi_size, dataSendReqs, MPI_STATUSES_IGNORE) != MPI_SUCCESS){
//					std::cout << "errors in communication (send)" << std::endl;
//				}
				double completeTime = this->myTimer_->stop();
				this->completeTimeMult_ += completeTime;
				double communicationTime = completeTime - computationTime;

//				if(sg::parallel::myGlobalMPIComm->getMyRank() == 1){
//					double sum_tmp = 0, sum_ref = 0;
//					double total_diff = 0;
//					int diffcounter = 0;
//					std::strstream indices;
//					for(int i = 0; i<this->numPatchedTrainingInstances_;i++){
//						if(fabs(temp[i] - temp_reference[i]) != 0){
//							sum_tmp += temp[i];
//							sum_ref += temp_reference[i];
//							total_diff += fabs(temp[i] - temp_reference[i]);
//							diffcounter++;
//							indices << i << " ";
//						}
//					}
//					if(total_diff != 0){
//						std::cout << "summs do not match: " << " " << sum_tmp << "!=" << sum_ref  << "  total diff: " << total_diff  << "; # of diffs: " << diffcounter<< std::endl << indices.str() << std::endl;
//						throw new sg::base::operation_exception("reuslts do not match!");
//					}
//					for(int i = 0; i<this->numPatchedTrainingInstances_;i++){

//						if(temp[i] - temp_reference[i] != 0.0){
//							std::cout << "results do not match: " << i << " " << temp[i] << "!=" << temp_reference[i]  << "diff" << temp[i] - temp_reference[i] << std::endl;
//							throw new sg::base::operation_exception("reuslts do not match!");
//						}
//					}
//				}

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

				/// @todo min/max mpireduce to see where the time is spent

				// patch result -> set additional entries zero
//				if (this->numTrainingInstances_ != temp.getSize())
//				{
//					for (size_t i = 0; i < (temp.getSize()-this->numTrainingInstances_); i++)
//					{
//						temp.set(temp.getSize()-(i+1), 0.0f);
//					}
//				}

				this->myTimer_->start();

			} //implicit OpenMP barrier here

			size_t myGridChunkStart = _mpi_grid_offsets_global[mpi_myrank];
			size_t myGridChunkEnd = myGridChunkStart + _mpi_grid_sizes_global[mpi_myrank];

			for(size_t thread_chunk = myGridChunkStart + thread_num; thread_chunk<myGridChunkEnd; thread_chunk+=thread_count){
				size_t start = _mpi_grid_offsets[thread_chunk];
				size_t end =  start + _mpi_grid_sizes[thread_chunk];

				MultTransType::multTranspose(level_, index_, dataset_, temp, result, start, end, 0, this->numPatchedTrainingInstances_);
				sg::parallel::myGlobalMPIComm->IsendToAll(&ptrResult[start], _mpi_grid_sizes[thread_chunk], start*2 + 3, &gridSendReqs[(thread_chunk - myGridChunkStart)*mpi_size]);
			}

#pragma omp barrier
		}
		this->computeTimeMultTrans_ += this->myTimer_->stop();

		MPI_Waitall(_chunkCountGrid, gridRecvReqs, MPI_STATUSES_IGNORE);
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

		MPI_Request gridRecvReqs[_chunkCountGrid]; //allocating a little more than necessary, otherwise complicated index computations needed
		int tags[_chunkCountGrid];
		for(int i = 0; i<_chunkCountGrid; i++){
			tags[i] = _mpi_grid_offsets[i]+1;
		}
		sg::parallel::myGlobalMPIComm->IrecvFromAll(ptrB, _mpi_grid_sizes_global, _mpi_grid_offsets_global, _mpi_grid_sizes, _mpi_grid_offsets, tags, _chunkCountGrid, gridRecvReqs);
		MPI_Request gridSendReqs[mpi_size * _mpi_grid_sizes_global[mpi_myrank]];

		int thread_count = 1;
		int thread_num = 0;
	#pragma omp parallel private(thread_count, thread_num)
		{
	#ifdef _OPENMP
			thread_count = omp_get_num_threads();
			thread_num = omp_get_thread_num();
	#endif
			size_t myGridChunkStart = _mpi_grid_offsets_global[mpi_myrank];
			size_t myGridChunkEnd = myGridChunkStart + _mpi_grid_sizes_global[mpi_myrank];
			for(size_t thread_chunk = myGridChunkStart + thread_num; thread_chunk<myGridChunkEnd; thread_chunk+=thread_count){
				size_t start = _mpi_grid_offsets[thread_chunk];
				size_t end =  start + _mpi_grid_sizes[thread_chunk];

				MultTransType::multTranspose(level_, index_, dataset_, myClasses, b, start, end, 0, this->numPatchedTrainingInstances_);
				//std::cout << "alltoall[" << thread_chunk << "]: from " << start << "; size" << _mpi_grid_sizes[thread_chunk] <<   std::endl;
				sg::parallel::myGlobalMPIComm->IsendToAll(&ptrB[start], _mpi_grid_sizes[thread_chunk], start+1, &gridSendReqs[(thread_chunk - myGridChunkStart)*mpi_size]);
				//std::cout << "after isendtoall " << thread_chunk << " --- " << _mpi_grid_offsets[thread_chunk] <<  std::endl;
			}
#pragma omp barrier
		}
		MPI_Status stats[_chunkCountGrid];
		MPI_Waitall(_chunkCountGrid, gridRecvReqs, stats);

		for(int i= 0; i<_chunkCountGrid; i++) {
			if(stats[i].MPI_ERROR != MPI_SUCCESS){
				std::cout << "["<< mpi_myrank <<"] status error at index "<< i << ": " << stats[i].MPI_ERROR << std::endl;
				throw new sg::base::operation_exception("DMSystemMatrixSPVectorizedIdentity : status error!");
			}
		}
	}

	virtual void rebuildLevelAndIndex(){
		delete this->level_;
		delete this->index_;

		this->level_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());
		this->index_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());

		m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

		int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
		size_t thread_count = 1;
#ifdef _OPENMP
#pragma omp parallel
		{
			thread_count = omp_get_num_threads();
		}
#endif
		size_t minChunkCount = mpi_size * thread_count; // every process should have at least thread_count chunks

		_chunkCountGrid = m_grid.getSize()/APPROXCHUNKSIZEGRID_ASYNC;
		_chunkCountGrid = std::max<size_t>(_chunkCountGrid, minChunkCount);


		delete[] _mpi_grid_sizes;
		delete[] _mpi_grid_offsets;

		// arrays for distribution settings
		_mpi_grid_sizes = new int[_chunkCountGrid];
		_mpi_grid_offsets = new int[_chunkCountGrid];

		sg::parallel::PartitioningTool::calcDistribution(m_grid.getSize(), _chunkCountGrid, _mpi_grid_sizes, _mpi_grid_offsets, 1);
		sg::parallel::PartitioningTool::calcDistribution(_chunkCountGrid, sg::parallel::myGlobalMPIComm->getNumRanks(), _mpi_grid_sizes_global, _mpi_grid_offsets_global, 1);

		if(sg::parallel::myGlobalMPIComm->getMyRank() == 0) {
			for(int i = 0; i<sg::parallel::myGlobalMPIComm->getNumRanks(); i++){
				//std::cout << "glob: " << _mpi_grid_offsets_global[i] << " (" << _mpi_grid_sizes_global[i] << ")" <<std::endl;
			}
			for (int i = 0; i<_chunkCountGrid; i++){
				//std::cout << "loc: " << _mpi_grid_offsets[i] << " (" << _mpi_grid_sizes[i] << ")" << std::endl;
			}
		}
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
};

}
}

#endif /* DMSYSTEMMATRIXVECTORIZEDIDENTITYASYNCMPI_HPP */

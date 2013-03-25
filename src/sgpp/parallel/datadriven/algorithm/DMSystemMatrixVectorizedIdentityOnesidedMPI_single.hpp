/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef DMSYSTEMMATRIXVECTORIZEDIDENTITYONESIDEDMPI_HPP
#define DMSYSTEMMATRIXVECTORIZEDIDENTITYONESIDEDMPI_HPP

#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityMPIBase.hpp"
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
template<typename KernelImplementation>
class DMSystemMatrixVectorizedIdentityOnesidedMPI : public sg::parallel::DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation::kernelType>
{
public:
	/**
	 * Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param trainData reference to sg::base::DataMatrix that contains the training data
	 * @param lambda the lambda, the regression parameter
	 * @param vecMode vectorization mode
	 */
	DMSystemMatrixVectorizedIdentityOnesidedMPI(
			sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
		: DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation::kernelType>(SparseGrid, trainData, lambda, vecMode),
		  _mpi_data_window(),
		  _mpi_data_window_buffer(NULL),
		  _mpi_grid_window(),
		  _mpi_grid_window_buffer(NULL)
	{
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

		_mpi_grid_sizes_global = new int[mpi_size];
		_mpi_grid_offsets_global = new int[mpi_size];
		_mpi_grid_sizes = NULL; // allocation in rebuildLevelAndIndex();
		_mpi_grid_offsets = NULL; // allocation in rebuildLevelAndIndex();
		rebuildLevelAndIndex();

		createAndInitDataBuffers();
		// mult: distribute calculations over dataset
		// multTranspose: distribute calculations over grid
	}

	void createAndInitGridBuffers(){
		if(_mpi_grid_window_buffer != NULL){
			delete _mpi_grid_window_buffer;
			MPI_Win_free(&_mpi_grid_window);
		}

		_mpi_grid_window_buffer = new sg::base::DataVector(this->storage_->size());

		MPI_Win_create(_mpi_grid_window_buffer->getPointer(), this->storage_->size()*sizeof(double),
					   sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &_mpi_grid_window);
		MPI_Win_fence(0, _mpi_grid_window);
	}

	void createAndInitDataBuffers(){
		// create MPI windows for databuffer
		if(_mpi_data_window_buffer != NULL){
			delete _mpi_data_window_buffer;
			MPI_Win_free(&_mpi_data_window);
		}
		_mpi_data_window_buffer = new sg::base::DataVector(this->numPatchedTrainingInstances_);
		MPI_Win_create(_mpi_data_window_buffer->getPointer(), this->numPatchedTrainingInstances_*sizeof(double),
					   sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &_mpi_data_window);
		MPI_Win_fence(0, _mpi_data_window);
	}

	/**
	 * Destructor
	 */
	virtual ~DMSystemMatrixVectorizedIdentityOnesidedMPI(){
		delete[] this->_mpi_grid_sizes;
		delete[] this->_mpi_grid_offsets;
		delete[] this->_mpi_data_sizes;
		delete[] this->_mpi_data_offsets;

		delete[] this->_mpi_grid_sizes_global;
		delete[] this->_mpi_grid_offsets_global;
		delete[] this->_mpi_data_sizes_global;
		delete[] this->_mpi_data_offsets_global;
		MPI_Win_free(&_mpi_data_window);
	}

	virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result){
		sg::base::DataVector temp(this->numPatchedTrainingInstances_);

		result.setAll(0.0);
		temp.setAll(0.0);
		_mpi_grid_window_buffer->setAll(0.0);
		_mpi_data_window_buffer->setAll(0.0);
		double* ptrResult = result.getPointer();
		double* ptrTemp = temp.getPointer();

		int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
		int mpi_myrank = sg::parallel::myGlobalMPIComm->getMyRank();

		MPI_Win_fence(0, _mpi_grid_window);
		MPI_Win_fence(0, _mpi_data_window);

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
				KernelImplementation::mult(
							this->level_,
							this->index_,
							this->mask_,
							this->offset_,
							this->dataset_,
							alpha,
							temp,
							0,
							alpha.getSize(),
							start,
							end);
				// patch result -> set additional entries zero
				// only done for processes that need this part of the temp data for multTrans
				for (size_t i = std::max<size_t>(this->numTrainingInstances_, start); i < end; i++)
				{
					temp.set(i, 0.0f);
				}
				myGlobalMPIComm->putToAll(&ptrTemp[start], start, _mpi_data_sizes[chunkIndex], _mpi_data_window);
			}

// this barrier is absolutely necessary, as it guarantees that all
// put-calls from all threads have been made. This is required as a MPI_Win_fence
// follows. If that fence is called before one or more threads have not finished their
// loop (which is very likely), the put won't be taken into account for the
// synchronization.
#pragma omp barrier
#pragma omp single
			{
				this->computeTimeMult_ += this->myTimer_->stop();

				MPI_Win_fence(0, _mpi_data_window);

				this->completeTimeMult_ += this->myTimer_->stop();
				this->myTimer_->start();

			} //implicit OpenMP barrier here

			size_t myGridChunkStart = _mpi_grid_offsets_global[mpi_myrank];
			size_t myGridChunkEnd = myGridChunkStart + _mpi_grid_sizes_global[mpi_myrank];
			sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
					myGridChunkStart, myGridChunkEnd, &threadStart, &threadEnd, 1);

			for(size_t thread_chunk = threadStart; thread_chunk<threadEnd; thread_chunk++){
				size_t start = _mpi_grid_offsets[thread_chunk];
				size_t end =  start + _mpi_grid_sizes[thread_chunk];

				KernelImplementation::multTranspose(
						this->level_,
						this->index_,
						this->mask_,
						this->offset_,
						this->dataset_,
						*_mpi_data_window_buffer,
						result,
						start,
						end,
						0,
						this->numPatchedTrainingInstances_);
				myGlobalMPIComm->putToAll(&ptrResult[start], start, _mpi_grid_sizes[thread_chunk], _mpi_grid_window);
			}
		}
		this->computeTimeMultTrans_ += this->myTimer_->stop();

		MPI_Win_fence(0, _mpi_grid_window);
		this->completeTimeMultTrans_ += this->myTimer_->stop();
		result.copyFrom(*_mpi_grid_window_buffer);
		result.axpy(static_cast<double>(this->numTrainingInstances_)*this->lambda_, alpha);
		if (mpi_myrank == 0) std::cout << "*";
	} //end mult

	virtual void generateb(sg::base::DataVector& classes, sg::base::DataVector& b){
		int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
		int mpi_myrank = sg::parallel::myGlobalMPIComm->getMyRank();

		_mpi_grid_window_buffer->setAll(0.0);
		MPI_Win_fence(0, _mpi_grid_window);

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
			size_t myGridChunkStart = _mpi_grid_offsets_global[mpi_myrank];
			size_t myGridChunkEnd = myGridChunkStart + _mpi_grid_sizes_global[mpi_myrank];
			size_t threadStart, threadEnd;
			sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
					myGridChunkStart, myGridChunkEnd, &threadStart, &threadEnd, 1);

			for(size_t thread_chunk = threadStart; thread_chunk<threadEnd; thread_chunk++){
				size_t start = _mpi_grid_offsets[thread_chunk];
				size_t end =  start + _mpi_grid_sizes[thread_chunk];

				KernelImplementation::multTranspose(
							this->level_,
							this->index_,
							this->mask_,
							this->offset_,
							this->dataset_,
							myClasses,
							b,
							start,
							end,
							0,
							this->numPatchedTrainingInstances_);
				myGlobalMPIComm->putToAll(&ptrB[start], start, _mpi_grid_sizes[thread_chunk], _mpi_grid_window);
			}
		}

		MPI_Win_fence(0, _mpi_grid_window);
		b.copyFrom(*_mpi_grid_window_buffer);
	}

	virtual void rebuildLevelAndIndex(){
		DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation::kernelType>::rebuildLevelAndIndex();

		if(_mpi_grid_sizes != NULL){
			delete[] _mpi_grid_sizes;
		}
		if(_mpi_grid_offsets != NULL){
			delete[] _mpi_grid_offsets;
		}

		/* calculate distribution of grid */
		calculateChunkCountGrid();
		_mpi_grid_sizes = new int[_chunkCountGrid];
		_mpi_grid_offsets = new int[_chunkCountGrid];
		sg::parallel::PartitioningTool::calcDistribution(this->storage_->size(), _chunkCountGrid, _mpi_grid_sizes, _mpi_grid_offsets, 1);
		sg::parallel::PartitioningTool::calcDistribution(_chunkCountGrid, sg::parallel::myGlobalMPIComm->getNumRanks(), _mpi_grid_sizes_global, _mpi_grid_offsets_global, 1);
		createAndInitGridBuffers();
	}

private:
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

	/// MPI windows
	sg::base::DataVector* _mpi_grid_window_buffer;
	sg::base::DataVector* _mpi_data_window_buffer;
	MPI_Win _mpi_grid_window;
	MPI_Win _mpi_data_window;

#define APPROXCHUNKSIZEGRID_ASYNC 20 // approximately how many blocks should be computed for the grid before sending
#define APPROXCHUNKSIZEDATA_ASYNC 150 // approximately how many blocks should be computed for the data before sending

	void calculateChunkCountGrid() {
		_chunkCountGrid = this->storage_->size()/APPROXCHUNKSIZEGRID_ASYNC;
		_chunkCountGrid = std::max<size_t>(_chunkCountGrid, getMinChunkCount());
		int chunkSize =
				static_cast<int>(
					std::ceil(
						static_cast<double>(this->storage_->size())/static_cast<double>(_chunkCountGrid)
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
#endif // DMSYSTEMMATRIXVECTORIZEDIDENTITYONESIDEDMPI_HPP

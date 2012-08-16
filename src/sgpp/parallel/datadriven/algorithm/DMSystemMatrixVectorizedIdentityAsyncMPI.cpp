/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/tools/MPI/SGppMPITools.hpp"
#include "base/exception/operation_exception.hpp"

#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/ChunkSizes.h"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityAsyncMPI.hpp"
#include "parallel/operation/ParallelOpFactory.hpp"
#include "parallel/tools/PartitioningTool.hpp"

namespace sg
{
namespace parallel
{

DMSystemMatrixVectorizedIdentityAsyncMPI::DMSystemMatrixVectorizedIdentityAsyncMPI(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
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

	size_t approximateChunkSizeData = 10; // approximately how many blocks should be computed before sending
	size_t blockCountData = this->numPatchedTrainingInstances_/sg::parallel::DMVectorizationPaddingAssistant::getVecWidthSP(this->vecMode_); // this process has (in total) blockCountData blocks of Data to process
	_chunkCountData = blockCountData/approximateChunkSizeData;

	// arrays for distribution settings
	_mpi_data_sizes = new int[_chunkCountData];
	_mpi_data_offsets = new int[_chunkCountData];
	_mpi_data_sizes_global = new int[mpi_size];
	_mpi_data_offsets_global = new int[mpi_size];

	calcDistribution(this->numPatchedTrainingInstances_, _chunkCountData, _mpi_data_sizes, _mpi_data_offsets, sg::parallel::DMVectorizationPaddingAssistant::getVecWidthSP(this->vecMode_));
	calcDistribution(_chunkCountData, mpi_size, _mpi_data_sizes_global, _mpi_data_offsets_global, 1);


	size_t approximateChunkSizeGrid = 100; // approximately how many blocks should be computed before sending
	_chunkCountGrid = m_grid.getStorage()->size()/approximateChunkSizeGrid;

	// arrays for distribution settings
	_mpi_grid_sizes = new int[_chunkCountGrid];
	_mpi_grid_offsets = new int[_chunkCountGrid];
	_mpi_grid_sizes_global = new int[mpi_size];
	_mpi_grid_offsets_global = new int[mpi_size];

	calcDistribution(m_grid.getStorage()->size(), _chunkCountGrid, _mpi_grid_sizes, _mpi_grid_offsets, 1);
	calcDistribution(m_grid.getStorage()->size(), mpi_size, _mpi_grid_sizes_global, _mpi_grid_offsets_global, 1);

	this->level_ = new sg::base::DataMatrix(m_grid.getStorage()->size(), m_grid.getStorage()->dim());
	this->index_ = new sg::base::DataMatrix(m_grid.getStorage()->size(), m_grid.getStorage()->dim());

	m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	// mult: distribute calculations over dataset
	// multTranspose: distribute calculations over grid
}

DMSystemMatrixVectorizedIdentityAsyncMPI::~DMSystemMatrixVectorizedIdentityAsyncMPI()
{
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

void DMSystemMatrixVectorizedIdentityAsyncMPI::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(this->numPatchedTrainingInstances_);
	result.setAll(0.0);
	double* ptrResult = result.getPointer();
	double* ptrTemp = temp.getPointer();

	int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
	int mpi_myrank = sg::parallel::myGlobalMPIComm->getMyRank();


	MPI_Request dataRecvReqs[mpi_size]; //allocating one more than necessary, otherwise complicated index computations needed
	MPI_Request gridRecvReqs[mpi_size]; //allocating one more than necessary, otherwise complicated index computations needed


	int rank;
	//posting temp reveices
	rank = 0;
	for(int i = 0; i<_chunkCountData; i++){
		// adjust rank to match chunk i
		if(i>=_mpi_data_offsets_global[rank] + _mpi_data_sizes_global[rank]){
			rank++;
		}
		// skip segments of this process, they are already there
		if(rank == mpi_myrank){
			i=_mpi_data_offsets_global[rank] + _mpi_data_sizes_global[rank]-1;
			// continue does the i++ (like after every iteration), so we are
			// at _mpi_data_offsets_global[rank] + _mpi_data_sizes_global[rank] at
			// the beginning of the next iteration, which means that we skipped rank mpi_myrank
			continue;
		}
		MPI_Irecv(&ptrTemp[_mpi_data_offsets[i]], _mpi_data_sizes[i], MPI_DOUBLE, rank, _mpi_data_offsets[i]*2+2, MPI_COMM_WORLD, &dataRecvReqs[rank]);
	}

	//posting grid receives
	rank = 0;
	for(int i = 0; i<_chunkCountGrid; i++){
		// adjust rank to match chunk i
		if(i>=_mpi_grid_offsets_global[rank] + _mpi_grid_sizes_global[rank]){
			rank++;
		}
		// skip segments of this process, they are already there
		if(rank == mpi_myrank){
			i=_mpi_grid_offsets_global[rank] + _mpi_grid_sizes_global[rank]-1;
			// continue does the i++ (like after every iteration), so we are
			// at _mpi_grid_offsets_global[rank] + _mpi_grid_sizes_global[rank] at
			// the beginning of the next iteration, which means that we skipped rank mpi_myrank
			continue;
		}
		MPI_Irecv(&ptrTemp[_mpi_grid_offsets[i]], _mpi_grid_sizes[i], MPI_DOUBLE, rank, _mpi_grid_offsets[i]*2+1, MPI_COMM_WORLD, &gridRecvReqs[rank]);
	}

	//mult
#ifdef _OPENMP
	#pragma omp parallel
	{
#endif
		size_t myDataChunkStart = _mpi_data_offsets_global[mpi_myrank];
		size_t myDataChunkEnd = myDataChunkStart + _mpi_data_sizes_global[mpi_myrank];
		for(size_t chunk = myDataChunkStart; chunk<myDataChunkEnd; chunk++){
			size_t start;
			size_t end;
			sg::parallel::PartitioningTool::getOpenMPLoopPartitionSegment(_mpi_data_offsets[chunk], _mpi_data_offsets[chunk] + _mpi_data_sizes[chunk], &start, &end, CHUNKDATAPOINTS_X86);

			if (start % CHUNKDATAPOINTS_X86 != 0 || end % CHUNKDATAPOINTS_X86 != 0)
			{
				std::cout << "start%CHUNKDATAPOINTS_X86: " << start%CHUNKDATAPOINTS_X86 << "; end%CHUNKDATAPOINTS_X86: " << end%CHUNKDATAPOINTS_X86 << std::endl;
				throw sg::base::operation_exception("processed vector segment must fit to CHUNKDATAPOINTS_X86!");
			}

			for(size_t c = start; c < end; c+=std::min<size_t>((size_t)CHUNKDATAPOINTS_X86, (end-c)))
			{
				//multV(this->level_, this->index_, this->dataset_, alpha, temp, c, end);
			}
			sg::parallel::myGlobalMPIComm->IsendToAll(&ptrTemp[_mpi_data_offsets[chunk]], _mpi_data_sizes[chunk], _mpi_data_offsets[chunk]*2+2);
		}

//	sg::parallel::myGlobalMPIComm->dataVectorAllToAll(temp, _mpi_grid_offsets, _mpi_grid_sizes);
#ifdef _OPENMP
	#pragma omp single
		{
#endif

	MPI_Waitall(mpi_size, dataRecvReqs, MPI_STATUSES_IGNORE);
	// patch result -> set additional entries zero
	if (this->numTrainingInstances_ != temp.getSize())
	{
		for (size_t i = 0; i < (temp.getSize()-this->numTrainingInstances_); i++)
		{
			temp.set(temp.getSize()-(i+1), 0.0f);
		}
	}
#ifdef _OPENMP
		} //implicit barrier here
#endif

	//multTrans
		size_t myGridChunkStart = _mpi_grid_offsets_global[mpi_myrank];
		size_t myGridChunkEnd = myDataChunkStart + _mpi_grid_sizes_global[mpi_myrank];
		for(size_t chunk = myGridChunkStart; chunk<myGridChunkEnd; chunk++){
			size_t start;
			size_t end;
			sg::parallel::PartitioningTool::getOpenMPLoopPartitionSegment(_mpi_grid_offsets[chunk], _mpi_grid_offsets[chunk] + _mpi_grid_offsets[chunk], &start, &end, 1);

			for(size_t k = start; k < end; k+=std::min<size_t>((size_t)CHUNKGRIDPOINTS_X86, (end-k)))
			{
				//multTransV(this->level_, this->index_, this->dataset_, temp, result, k, end);
			}
			sg::parallel::myGlobalMPIComm->IsendToAll(&ptrResult[_mpi_grid_offsets[chunk]], _mpi_grid_sizes[chunk], _mpi_grid_offsets[chunk]*2+1);
		}

#ifdef _OPENMP
	}
#endif
	MPI_Waitall(mpi_size, gridRecvReqs, MPI_STATUSES_IGNORE);

	result.axpy(static_cast<double>(this->numTrainingInstances_)*this->lambda_, alpha);
}

void DMSystemMatrixVectorizedIdentityAsyncMPI::generateb(sg::base::DataVector& classes, sg::base::DataVector& b)
{
	sg::base::DataVector myClasses(classes);

	// Apply padding
	if (this->numPatchedTrainingInstances_ != myClasses.getSize())
	{
		myClasses.resizeZero(this->numPatchedTrainingInstances_);
	}

	multTransposeVec(myClasses, b);

	debugMPI(sg::parallel::myGlobalMPIComm, "end generate b");
}

void DMSystemMatrixVectorizedIdentityAsyncMPI::rebuildLevelAndIndex()
{
	delete this->level_;
	delete this->index_;

	this->level_ = new sg::base::DataMatrix(m_grid.getStorage()->size(), m_grid.getStorage()->dim());
	this->index_ = new sg::base::DataMatrix(m_grid.getStorage()->size(), m_grid.getStorage()->dim());

	m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	size_t approximateChunkSizeGrid = 100; // approximately how many blocks should be computed before sending
	_chunkCountGrid = m_grid.getStorage()->size()/approximateChunkSizeGrid;

	delete[] _mpi_grid_sizes;
	delete[] _mpi_grid_offsets;

	// arrays for distribution settings
	_mpi_grid_sizes = new int[_chunkCountGrid];
	_mpi_grid_offsets = new int[_chunkCountGrid];

	calcDistribution(m_grid.getStorage()->size(), _chunkCountGrid, _mpi_grid_sizes, _mpi_grid_offsets, 1);
	calcDistribution(m_grid.getStorage()->size(), sg::parallel::myGlobalMPIComm->getNumRanks(), _mpi_grid_sizes_global, _mpi_grid_offsets_global, 1);

}

void DMSystemMatrixVectorizedIdentityAsyncMPI::multVec(base::DataVector &alpha, base::DataVector &result)
{
	this->myTimer_->start();

	//this->computeTimeMult_ += this->B_->multVectorized(alpha, result);

	debugMPI(sg::parallel::myGlobalMPIComm, "_mpi_data_offsets[" << sg::parallel::myGlobalMPIComm->getMyRank() << "] = " << _mpi_data_offsets[sg::parallel::myGlobalMPIComm->getMyRank()] << std::endl);
	debugMPI(sg::parallel::myGlobalMPIComm, "_mpi_data_sizes[" << sg::parallel::myGlobalMPIComm->getMyRank() << "] = " << _mpi_data_sizes[sg::parallel::myGlobalMPIComm->getMyRank()] << std::endl);

	debugMPI(sg::parallel::myGlobalMPIComm, "result.size() = " << result.getSize());

	sg::parallel::myGlobalMPIComm->dataVectorAllToAll(result, _mpi_data_offsets, _mpi_data_sizes);

	this->completeTimeMult_ += this->myTimer_->stop();
}

void DMSystemMatrixVectorizedIdentityAsyncMPI::multTransposeVec(base::DataVector &source, base::DataVector &result)
{
	this->myTimer_->start();
	//this->computeTimeMultTrans_ += this->B_->multTransposeVectorized(source, result);

	sg::parallel::myGlobalMPIComm->dataVectorAllToAll(result, _mpi_grid_offsets, _mpi_grid_sizes);

	this->completeTimeMultTrans_ += this->myTimer_->stop();
}

void DMSystemMatrixVectorizedIdentityAsyncMPI::calcDistribution(int totalSize, int numChunks, int *sizes, int *offsets, size_t blocksize)
{
	for(int chunkID = 0; chunkID < numChunks; ++chunkID){
		size_t size;
		size_t offset;
		sg::parallel::PartitioningTool::getPartitionSegment(totalSize, numChunks, chunkID, &size, &offset, blocksize);
		sizes[chunkID] = size;
		offsets[chunkID] = offset;
	}
}

}
}

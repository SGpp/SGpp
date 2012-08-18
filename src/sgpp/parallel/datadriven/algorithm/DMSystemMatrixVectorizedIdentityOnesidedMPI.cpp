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
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinearMult.h"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinearMultTranspose.h"
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

	int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
	int mpi_myrank = sg::parallel::myGlobalMPIComm->getMyRank();

	sg::parallel::X86SimdLinearMult mult = sg::parallel::X86SimdLinearMult(level_, index_, dataset_, alpha, temp);
	double* ptrTemp = temp.getPointer();
	MPI_Win tempWin;
	MPI_Win_create(ptrTemp, sizeof(double) * temp.getSize(), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &tempWin);

	sg::parallel::X86SimdLinearMultTranspose multTranspose = sg::parallel::X86SimdLinearMultTranspose(level_, index_, dataset_, temp, result);
	double* ptrResult = result.getPointer();
	MPI_Win resultWin;
	MPI_Win_create(ptrResult, sizeof(double) * result.getSize(), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &resultWin);
	this->myTimer_->start();

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

			mult(start, end);
			sg::parallel::myGlobalMPIComm->putToAll(temp, start, end-start, tempWin);
		}
#ifdef _OPENMP
	#pragma omp single
		{
#endif

	MPI_Waitall(mpi_size, dataRecvReqs, MPI_STATUSES_IGNORE);
	this->completeTimeMult_ += this->myTimer_->stop();

	// patch result -> set additional entries zero
	if (this->numTrainingInstances_ != temp.getSize())
	{
		for (size_t i = 0; i < (temp.getSize()-this->numTrainingInstances_); i++)
		{
			temp.set(temp.getSize()-(i+1), 0.0f);
		}
	}

	this->myTimer_->start();

#ifdef _OPENMP
		} //implicit OpenMP barrier here
#endif

		multTransposeVec(multTranspose);

#ifdef _OPENMP
	}
#endif
	MPI_Waitall(mpi_size, gridRecvReqs, MPI_STATUSES_IGNORE);
	this->completeTimeMultTrans_ += this->myTimer_->stop();

	result.axpy(static_cast<double>(this->numTrainingInstances_)*this->lambda_, alpha);
}

void DMSystemMatrixVectorizedIdentityAsyncMPI::generateb(sg::base::DataVector& classes, sg::base::DataVector& b)
{
	int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();

	sg::base::DataVector myClasses(classes);
	// Apply padding
	if (this->numPatchedTrainingInstances_ != myClasses.getSize())
	{
		myClasses.resizeZero(this->numPatchedTrainingInstances_);
	}

	MPI_Request gridRecvReqs[mpi_size]; //allocating one more than necessary, otherwise complicated index computations needed
	sg::parallel::myGlobalMPIComm->IrecvFromAll(b.getPointer(), _mpi_grid_sizes_global, _mpi_grid_offsets_global, _mpi_grid_sizes, _mpi_grid_offsets, _chunkCountGrid, gridRecvReqs);

	sg::parallel::X86SimdLinearMultTranspose multTranspose = sg::parallel::X86SimdLinearMultTranspose(level_, index_, dataset_, myClasses, b);

	std::cout << "before openmp" << std::endl;

#ifdef _OPENMP
	#pragma omp parallel
	{
#endif
		std::cout << "in openmp" << std::endl;

		multTransposeVec(multTranspose);
		std::cout << "in openmp 2" << std::endl;

#ifdef _OPENMP
	}
#endif

	std::cout << "after openmp" << std::endl;

	MPI_Waitall(mpi_size, gridRecvReqs, MPI_STATUSES_IGNORE);

	std::cout << "after waitall" << std::endl;
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

void DMSystemMatrixVectorizedIdentityAsyncMPI::multVec(sg::parallel::X86SimdLinearMult &mult)
{

}

void DMSystemMatrixVectorizedIdentityAsyncMPI::multTransposeVec(sg::parallel::X86SimdLinearMultTranspose &multTranspose)
{
	int mpi_myrank = sg::parallel::myGlobalMPIComm->getMyRank();

	double* ptrResult = multTranspose._result.getPointer();

	size_t myGridChunkStart = _mpi_grid_offsets_global[mpi_myrank];
	size_t myGridChunkEnd = myGridChunkStart + _mpi_grid_sizes_global[mpi_myrank];
	for(size_t chunk = myGridChunkStart; chunk<myGridChunkEnd; chunk++){
		size_t start;
		size_t end;
		sg::parallel::PartitioningTool::getOpenMPLoopPartitionSegment(_mpi_grid_offsets[chunk], _mpi_grid_offsets[chunk] + _mpi_grid_offsets[chunk], &start, &end, 1);

		multTranspose(start, end);
		std::cout << "before isendtoall" << std::endl;

		sg::parallel::myGlobalMPIComm->IsendToAll(&ptrResult[start], end-start);
		std::cout << "after isendtoall" << _mpi_grid_offsets[chunk] <<  std::endl;
	}
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

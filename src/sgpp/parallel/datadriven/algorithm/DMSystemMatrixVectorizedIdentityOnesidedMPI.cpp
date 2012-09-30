/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityOnesidedMPI.hpp"
#include "base/exception/operation_exception.hpp"

#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinearMult.h"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinearMultTranspose.h"
#include "parallel/operation/ParallelOpFactory.hpp"
#include "parallel/tools/PartitioningTool.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

#define APPROXCHUNKSIZEGRID 50 // approximately how many blocks should be computed for the grid before sending
#define APPROXCHUNKSIZEDATA 50 // approximately how many blocks should be computed for the data before sending

namespace sg
{
namespace parallel
{

DMSystemMatrixVectorizedIdentityOneSidedMPI::DMSystemMatrixVectorizedIdentityOneSidedMPI(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
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
	int thread_count = 1;
#ifdef _OPENMP
#pragma omp parallel
	{
		thread_count = omp_get_num_threads();
	}
#endif

	size_t blockCountData = this->numPatchedTrainingInstances_/sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_); // this process has (in total) blockCountData blocks of Data to process
	_chunkCountData = blockCountData/APPROXCHUNKSIZEDATA;

	// arrays for distribution settings
	_mpi_data_sizes = new int[_chunkCountData];
	_mpi_data_offsets = new int[_chunkCountData];
	_mpi_data_sizes_global = new int[mpi_size];
	_mpi_data_offsets_global = new int[mpi_size];

	calcDistribution(this->numPatchedTrainingInstances_, _chunkCountData, _mpi_data_sizes, _mpi_data_offsets, sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));
	calcDistribution(_chunkCountData, mpi_size, _mpi_data_sizes_global, _mpi_data_offsets_global, 1);

	if(sg::parallel::myGlobalMPIComm->getMyRank() == 0){
		for(int i = 0; i<125; i++){
			std::cout << "data_sizes[" << i << "]: " << _mpi_data_sizes[i] << std::endl;
			std::cout << "data_offsets[" << i << "]: " << _mpi_data_offsets[i] << std::endl;
		}
		for(int i = 0; i<mpi_size; i++){
			std::cout << "data_sizes_g[" << i << "]: " << _mpi_data_sizes_global[i] << std::endl;
			std::cout << "data_offsets_g[" << i << "]: " << _mpi_data_offsets_global[i] << std::endl;
		}
	}
	_chunkCountGrid = m_grid.getSize()/APPROXCHUNKSIZEGRID;

	// arrays for distribution settings
	_mpi_grid_sizes = new int[_chunkCountGrid];
	_mpi_grid_offsets = new int[_chunkCountGrid];
	_mpi_grid_sizes_global = new int[mpi_size];
	_mpi_grid_offsets_global = new int[mpi_size];

	calcDistribution(m_grid.getSize(), _chunkCountGrid, _mpi_grid_sizes, _mpi_grid_offsets, 1);
	calcDistribution(_chunkCountGrid, mpi_size, _mpi_grid_sizes_global, _mpi_grid_offsets_global, 1);


	this->level_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());
	this->index_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());

	m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	// create MPI windows
	_mpi_data_window_buffer = new sg::base::DataVector(this->numPatchedTrainingInstances_);
	_mpi_grid_window_buffer = new sg::base::DataVector(m_grid.getSize());

	double* ptrData = _mpi_data_window_buffer->getPointer();
	double* ptrGrid = _mpi_grid_window_buffer->getPointer();

	_mpi_data_window = new MPI_Win[mpi_size];
	_mpi_grid_window = new MPI_Win[mpi_size];
	for(int rank = 0; rank<mpi_size; rank++) {
		// data
		int size = 0;
		for(int i = 0; i<_mpi_data_sizes_global[rank]; i++){
			size += _mpi_data_sizes[_mpi_data_offsets_global[rank] + i];
		}
		int offset = _mpi_data_offsets[_mpi_data_offsets_global[rank]];
		MPI_Win_create(&ptrData[offset], size*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &_mpi_data_window[rank]);
		if(sg::parallel::myGlobalMPIComm->getMyRank() == 0)
			std::cout << "win init data [" << rank << "] " <<  offset << " - " << offset + size - 1 << " = " << size << std::endl;
		MPI_Win_fence(0, _mpi_data_window[rank]);

		// grid
		size = 0;
		for(int i = 0; i<_mpi_grid_sizes_global[rank]; i++){
			size += _mpi_grid_sizes[_mpi_grid_offsets_global[rank] + i];
		}
		offset = _mpi_grid_offsets[_mpi_grid_offsets_global[rank]];
		MPI_Win_create(&ptrGrid[offset], size*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &_mpi_grid_window[rank]);
		MPI_Win_fence(0, _mpi_grid_window[rank]);

		if(sg::parallel::myGlobalMPIComm->getMyRank() == 0)
			std::cout << "win init grid [" << rank << "] " <<  offset << " - " << offset + size - 1 << " = " << size << std::endl;
	}

	if(sg::parallel::myGlobalMPIComm->getMyRank() == 0 || true){
		for(int rank = 0; rank<mpi_size; rank++) {
			int base_addr, size, disp_unit, flag;
			MPI_Win *w = &_mpi_data_window[rank];
			MPI_Win_get_attr(*w, MPI_WIN_BASE, &base_addr, &flag);
			MPI_Win_get_attr(*w, MPI_WIN_SIZE, &size, &flag);
			MPI_Win_get_attr(*w, MPI_WIN_DISP_UNIT, &disp_unit, &flag);
			size_t offset = ((size_t)base_addr - (size_t)((void*)ptrData));
			size_t real_offset = offset/sizeof(double);
			size_t real_size = size/sizeof(double);
			std::cout << "[" << rank << "][mpi_data_window] offset:" <<  offset << " - realOffset:" <<  real_offset <<" - " << real_size << " . " << disp_unit << std::endl;
		}

		std::cout << "bla" << std::endl;
		for(int rank = 0; rank<mpi_size; rank++) {
			int base_addr, size, disp_unit, flag;
			MPI_Win *w = &_mpi_grid_window[rank];
			MPI_Win_get_attr(*w, MPI_WIN_BASE, &base_addr, &flag);
			if (!flag) {
				fprintf( stderr, "Attribute for key MPI_WIN_BASE not set\n" );
			}
			MPI_Win_get_attr(*w, MPI_WIN_SIZE, &size, &flag);
			if (!flag) {
				fprintf( stderr, "Attribute for key MPI_WIN_SIZE not set\n" );
			}
			MPI_Win_get_attr(*w, MPI_WIN_DISP_UNIT, &disp_unit, &flag);
			if (!flag) {
				fprintf( stderr, "Attribute for key MPI_WIN_DISP_UNIT not set\n" );
			}
			size_t offset = ((size_t)base_addr - (size_t)((void*)ptrGrid));
			size_t real_offset = offset/sizeof(double);
			disp_unit = size - disp_unit;
			size_t real_size = ((size_t)size)/(size_t)disp_unit;//sizeof(double);
			std::cout << "[" << rank << "][mpi_grid_window]" <<  offset << " - " <<  real_offset << " - real size: " << real_size << " - size: " << size  << " - disp_unit: " << disp_unit  << " - base: " << base_addr << std::endl;
		}

		std::cout << "bla" << std::endl;
		std::cout.flush();
	}
	MPI_Barrier(MPI_COMM_WORLD);



	// mult: distribute calculations over dataset
	// multTranspose: distribute calculations over grid
}

DMSystemMatrixVectorizedIdentityOneSidedMPI::~DMSystemMatrixVectorizedIdentityOneSidedMPI()
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

	int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
	for (int rank = 0; rank<mpi_size; rank++) {
		MPI_Win_free(&_mpi_grid_window[rank]);
		MPI_Win_free(&_mpi_data_window[rank]);
	}
}

void DMSystemMatrixVectorizedIdentityOneSidedMPI::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
	int mpi_myrank = sg::parallel::myGlobalMPIComm->getMyRank();
	for (int rank = 0; rank < mpi_size; rank++){
		int assert = MPI_MODE_NOPRECEDE;
		if(rank == mpi_myrank) {
			assert |= MPI_MODE_NOPUT;
		}
		MPI_Win_fence(0, _mpi_data_window[rank]);
		MPI_Win_fence(0, _mpi_grid_window[rank]);
	}
	result.setAll(0.0);
	_mpi_data_window_buffer->setAll(0.0);
	_mpi_grid_window_buffer->setAll(0.0);

	this->myTimer_->start();

	int thread_count = 1;
	int thread_num = 0;
#ifdef _OPENMP
	#pragma omp parallel private (thread_num, thread_count)
	{
		thread_count = omp_get_num_threads();
		thread_num = omp_get_thread_num();
#endif

		size_t myDataChunkStart = _mpi_data_offsets_global[mpi_myrank];
		size_t myDataChunkEnd = myDataChunkStart + _mpi_data_sizes_global[mpi_myrank];
		size_t procDataChunkStart = _mpi_data_offsets[myDataChunkStart];
		for(size_t thread_chunk = myDataChunkStart + thread_num; thread_chunk<myDataChunkEnd; thread_chunk+=thread_count){
			size_t start = _mpi_data_offsets[thread_chunk];
			size_t end = start + _mpi_data_sizes[thread_chunk];

			if (start % sg::parallel::X86SimdLinearMult::getChunkDataPoints() != 0 || end % sg::parallel::X86SimdLinearMult::getChunkDataPoints() != 0 || end > this->numPatchedTrainingInstances_)
			{
				std::cout << "start%CHUNKDATAPOINTS_X86: " << start%sg::parallel::X86SimdLinearMult::getChunkDataPoints() << "; end%CHUNKDATAPOINTS_X86: " << end%sg::parallel::X86SimdLinearMult::getChunkDataPoints() << std::endl;
				std::cout << end << " > numpatchedtraininstances ("<< this->numPatchedTrainingInstances_ <<")" << std::endl;
				throw sg::base::operation_exception("processed vector segment must fit to CHUNKDATAPOINTS_X86!");
			}
			sg::parallel::X86SimdLinearMult::mult(level_, index_, dataset_, alpha, *_mpi_data_window_buffer, start, end);
			sg::parallel::myGlobalMPIComm->putToAllInplace(_mpi_data_window[mpi_myrank], start - procDataChunkStart, end-start);
		}

#ifdef _OPENMP
	#pragma omp single
		{
#endif
	//std::cout << "[" <<  mpi_myrank << "] start mult " << std::endl;
	this->computeTimeMult_ += this->myTimer_->stop();
	for (int rank = 0; rank < mpi_size; rank++){
		int assert = MPI_MODE_NOSUCCEED;
		if(rank != mpi_myrank){
			assert |= MPI_MODE_NOSTORE;
		}
		MPI_Win_fence(assert, _mpi_data_window[rank]);
	}
	this->completeTimeMult_ += this->myTimer_->stop();

	// patch result -> set additional entries zero
	if (this->numTrainingInstances_ != _mpi_data_window_buffer->getSize())
	{
		for (size_t i = 0; i < (_mpi_data_window_buffer->getSize()-this->numTrainingInstances_); i++)
		{
			_mpi_data_window_buffer->set(_mpi_data_window_buffer->getSize()-(i+1), 0.0);
		}
	}

	this->myTimer_->start();

#ifdef _OPENMP
		} //implicit OpenMP barrier here
#endif

		size_t myGridChunkStart = _mpi_grid_offsets_global[mpi_myrank];
		size_t myGridChunkEnd = myGridChunkStart + _mpi_grid_sizes_global[mpi_myrank];
		size_t procGridChunkStart = _mpi_grid_offsets[myGridChunkStart];
		for(size_t thread_chunk = myGridChunkStart + thread_num; thread_chunk<myGridChunkEnd; thread_chunk+=thread_count){
			size_t start = _mpi_grid_offsets[thread_chunk];
			size_t end =  start + _mpi_grid_sizes[thread_chunk];

			sg::parallel::X86SimdLinearMultTranspose::multTranspose(level_, index_, dataset_, *_mpi_data_window_buffer, *_mpi_grid_window_buffer, start, end, 0, this->numPatchedTrainingInstances_);
			sg::parallel::myGlobalMPIComm->putToAllInplace(_mpi_grid_window[mpi_myrank], start - procGridChunkStart, end-start);
		}


#ifdef _OPENMP
	}
#endif
	this->computeTimeMultTrans_ += this->myTimer_->stop();

	for (int rank = 0; rank < mpi_size; rank++){
		int assert = MPI_MODE_NOSUCCEED;
		if(rank != mpi_myrank){
			assert |= MPI_MODE_NOSTORE;
		}
		MPI_Win_fence(assert, _mpi_grid_window[rank]);
	}

	this->completeTimeMultTrans_ += this->myTimer_->stop();
	result.copyFrom(*_mpi_grid_window_buffer);

	result.axpy(static_cast<double>(this->numTrainingInstances_)*this->lambda_, alpha);
	if (mpi_myrank == 0) std::cout << "*";
}

void DMSystemMatrixVectorizedIdentityOneSidedMPI::generateb(sg::base::DataVector& classes, sg::base::DataVector& b)
{
	int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();
	int mpi_myrank = sg::parallel::myGlobalMPIComm->getMyRank();

	// reuse window buffer, it already has the correct size
	_mpi_data_window_buffer->setAll(0.0);
	_mpi_data_window_buffer->copyFrom(classes);

	_mpi_grid_window_buffer->setAll(0.0);
	b.setAll(0.0);

	for (int rank = 0; rank < mpi_size; rank++){
		MPI_Win_fence(0, _mpi_grid_window[rank]);
	}

	int thread_count = 1;
	int thread_num = 0;
#ifdef _OPENMP
	#pragma omp parallel private (thread_num, thread_count)
	{
		thread_count = omp_get_num_threads();
		thread_num = omp_get_thread_num();
#endif

		size_t myGridChunkStart = _mpi_grid_offsets_global[mpi_myrank];
		size_t myGridChunkEnd = myGridChunkStart + _mpi_grid_sizes_global[mpi_myrank];
		size_t procChunkStart = _mpi_grid_offsets[myGridChunkStart];
		for(size_t thread_chunk = myGridChunkStart + thread_num; thread_chunk<myGridChunkEnd; thread_chunk+=thread_count){
			size_t start = _mpi_grid_offsets[thread_chunk];
			size_t end =  start + _mpi_grid_sizes[thread_chunk];

			sg::parallel::X86SimdLinearMultTranspose::multTranspose(level_, index_, dataset_, *_mpi_data_window_buffer, *_mpi_grid_window_buffer, start, end, 0, this->numPatchedTrainingInstances_);
			sg::parallel::myGlobalMPIComm->putToAllInplace(_mpi_grid_window[mpi_myrank], start - procChunkStart, end - start);
		}
#ifdef _OPENMP
	}
#endif
	for (int rank = 0; rank < mpi_size; rank++){
		MPI_Win_fence(0, _mpi_grid_window[rank]);
	}
	b.copyFrom(*_mpi_grid_window_buffer);

	MPI_Barrier(MPI_COMM_WORLD);
}

void DMSystemMatrixVectorizedIdentityOneSidedMPI::rebuildLevelAndIndex()
{
	delete this->level_;
	delete this->index_;

	this->level_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());
	this->index_ = new sg::base::DataMatrix(m_grid.getSize(), m_grid.getStorage()->dim());

	m_grid.getStorage()->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

	_chunkCountGrid = m_grid.getSize()/APPROXCHUNKSIZEGRID;

	delete[] _mpi_grid_sizes;
	delete[] _mpi_grid_offsets;

	// arrays for distribution settings
	_mpi_grid_sizes = new int[_chunkCountGrid];
	_mpi_grid_offsets = new int[_chunkCountGrid];

	int mpi_size = sg::parallel::myGlobalMPIComm->getNumRanks();

	calcDistribution(m_grid.getSize(), _chunkCountGrid, _mpi_grid_sizes, _mpi_grid_offsets, 1);
	calcDistribution(_chunkCountGrid, mpi_size, _mpi_grid_sizes_global, _mpi_grid_offsets_global, 1);


	_mpi_grid_window_buffer = new sg::base::DataVector(m_grid.getSize());
	double* ptrGrid = _mpi_grid_window_buffer->getPointer();

	for (int rank = 0; rank<mpi_size; rank++) {
		MPI_Win_free(&_mpi_grid_window[rank]);
	}

	_mpi_grid_window = new MPI_Win[mpi_size];
	for(int rank = 0; rank<mpi_size; rank++) {
		// grid
		int size = 0;
		for(int i = 0; i<_mpi_grid_sizes_global[rank]; i++){
			size += _mpi_grid_sizes[_mpi_grid_offsets_global[rank] + i];
		}
		int offset = _mpi_grid_offsets[_mpi_grid_offsets_global[rank]];
		MPI_Win_create(&ptrGrid[offset], size*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &_mpi_grid_window[rank]);
		MPI_Win_fence(0, _mpi_grid_window[rank]);

		if(sg::parallel::myGlobalMPIComm->getMyRank() == 0)
			std::cout << "win init grid [" << rank << "] " <<  offset << " - " << offset + size - 1 << " = " << size << std::endl;
	}
}


void DMSystemMatrixVectorizedIdentityOneSidedMPI::calcDistribution(int totalSize, int numChunks, int *sizes, int *offsets, size_t blocksize)
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
